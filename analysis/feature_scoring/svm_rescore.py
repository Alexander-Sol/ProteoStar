#!/usr/bin/env python
"""Semi-supervised SVM feature rescoring (feature-level Percolator).

Learns a multivariate quality score that separates real MS1 features from junk, using
target-decoy training:
  - NEGATIVES: features detected under decoy models (spacing / weird-averagine combs) -- these are
    detections the pipeline makes by chance under a wrong isotope model, i.e. labeled junk.
  - POSITIVES: a confident subset of TARGET features. Seeded by apex envelope fit (Decon Score),
    then re-selected each iteration from the current model's own out-of-fold score at a decoy-FDR
    threshold -- so the positive set evolves (the "semi-supervised" part).

3-fold cross-validation (train on 2 folds, score the held-out fold) breaks the self-training
overfit loop: a target is always scored by a model that never saw it, and positive re-selection
uses those out-of-fold scores. This mirrors Percolator (Kall 2007).

Win condition is NOT decoy AUC (leakage from seeding on a fit input can inflate it) but PSM
recall at matched feature count: does ranking targets by the SVM score retain recall better than
ranking by Summed Intensity? The eval prints both curves side by side.

Usage:
  svm_rescore.py --target T.tsv --decoys D1.tsv,D2.tsv[,...]
                 [--gt GT.tsv --rt 0.5] [--iters 10] [--seed-frac 0.5] [--q 0.05]
"""
import argparse, csv, bisect
import numpy as np
from sklearn.svm import LinearSVC
from sklearn.preprocessing import StandardScaler

# Quality/shape features only -- deliberately NO raw mass/RT/mz (the SVM must not exploit where in
# mass/RT space the decoy pass happened to land; it must learn feature *quality*).
FEATURE_NAMES = [
    "log_intensity", "num_charge_states", "primary_charge", "cross_charge_support",
    "num_members", "decon_score", "min_decon_score", "max_num_isotopes",
    "ppm_spread", "ppm_missing", "rt_extent",
]
DECON_IDX = FEATURE_NAMES.index("decon_score")
INTEN_IDX = FEATURE_NAMES.index("log_intensity")
MASS_PPM = 20.0


def load(path):
    """-> (X [n,d] float32, meta list of (mass, rt_apex, set(charges)))."""
    X, meta = [], []
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                si = float(row["Summed Intensity"])
                ncs = float(row["Num Charge States"])
                pc = float(row["Primary Charge"])
                ccs = float(row["Cross-Charge Support"])
                nm = float(row["Num Members"])
                ds = float(row["Decon Score"])
                mds = float(row["Min Decon Score"])
                mi = float(row["Max Num Isotopes"])
                ppm = float(row["PPM Spread"])
                rts = float(row["RT Start"]); rte = float(row["RT End"])
                mass = float(row["Monoisotopic Mass"]); rtap = float(row["RT Apex"])
            except (KeyError, ValueError):
                continue
            ppm_missing = 1.0 if ppm >= 999 else 0.0
            ppm_val = 15.0 if ppm >= 999 else ppm      # cap the single-isotope sentinel
            X.append([np.log1p(si), ncs, pc, ccs, nm, ds, mds, mi, ppm_val, ppm_missing, rte - rts])
            charges = set()
            for c in row.get("Charge States", "").split(";"):
                if c.strip():
                    charges.add(int(c))
            meta.append((mass, rtap, charges))
    return np.asarray(X, dtype=np.float64), meta


def recall_at(mask, meta, refs, rt_delta):
    """PSM recall using only target features where mask is True."""
    feats = sorted((meta[i] for i in np.nonzero(mask)[0]), key=lambda m: m[0])
    masses = [f[0] for f in feats]
    matched = 0
    for (rmass, rrt, _rz) in refs:
        tol = rmass * MASS_PPM / 1e6
        lo = bisect.bisect_left(masses, rmass - tol)
        hi = bisect.bisect_right(masses, rmass + tol)
        for i in range(lo, hi):
            if abs(feats[i][1] - rrt) <= rt_delta:
                matched += 1
                break
    return matched


def load_refs(path):
    refs = []
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                refs.append((float(row["mono_mass"]), float(row["rt"]), int(row["charge"])))
            except (ValueError, KeyError):
                continue
    return refs


def fit_score(Xtr, ytr, Xall):
    """Standardize on the training rows, fit a class-balanced linear SVM, return decision scores
    for every row in Xall."""
    sc = StandardScaler().fit(Xtr)
    clf = LinearSVC(C=1.0, class_weight="balanced", dual="auto", max_iter=20000)
    clf.fit(sc.transform(Xtr), ytr)
    return clf.decision_function(sc.transform(Xall))


def qvalues(scores_t, scores_d):
    """Decoy-estimated q-value for each target, treating pooled decoys as the null. For threshold s,
    FDR(s) = (frac of decoys >= s) * n_target / (# targets >= s); q = running min from high scores."""
    nt, nd = len(scores_t), max(len(scores_d), 1)
    order = np.argsort(-scores_t)
    ds_sorted = np.sort(scores_d)
    q = np.ones(nt)
    running = 1.0
    for rank, idx in enumerate(order, start=1):
        s = scores_t[idx]
        n_dec_ge = nd - bisect.bisect_left(ds_sorted, s)
        fdr = (n_dec_ge / nd) * nt / rank
        running = min(running, fdr) if rank == 1 else running
        q[idx] = fdr
    # enforce monotone q (running min from best score downward)
    qm = np.ones(nt); m = 1.0
    for idx in order:
        m = min(m, q[idx]); qm[idx] = m
    return qm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--decoys", required=True, help="comma-separated decoy TSVs")
    ap.add_argument("--gt"); ap.add_argument("--rt", type=float, default=0.5)
    ap.add_argument("--iters", type=int, default=10)
    ap.add_argument("--seed-frac", type=float, default=0.5, help="top fraction by Decon Score as initial positives")
    ap.add_argument("--q", type=float, default=0.05, help="decoy-FDR to admit positives each iteration")
    ap.add_argument("--folds", type=int, default=3)
    a = ap.parse_args()

    Xt, meta_t = load(a.target)
    Xd_list = [load(p)[0] for p in a.decoys.split(",")]
    Xd = np.vstack(Xd_list)
    nt, nd = len(Xt), len(Xd)
    print(f"target features: {nt}   decoy features: {nd}   ({', '.join(a.decoys.split(','))})")

    rng = np.random.default_rng(0)
    fold_t = rng.integers(0, a.folds, nt)
    fold_d = rng.integers(0, a.folds, nd)

    # Seed positives: top seed-frac of targets by Decon Score.
    thr = np.quantile(Xt[:, DECON_IDX], 1.0 - a.seed_frac)
    pos = Xt[:, DECON_IDX] >= thr
    print(f"seed positives (Decon Score >= p{100*(1-a.seed_frac):.0f} = {thr:.3f}): {pos.sum()}")

    oof_t = np.zeros(nt)   # out-of-fold target scores
    for it in range(a.iters):
        oof_t[:] = 0.0
        oof_d = np.zeros(nd)
        for f in range(a.folds):
            tr_t = (fold_t != f) & pos          # positive targets NOT in this fold
            tr_d = fold_d != f                  # decoys NOT in this fold
            Xtr = np.vstack([Xt[tr_t], Xd[tr_d]])
            ytr = np.concatenate([np.ones(tr_t.sum()), np.zeros(tr_d.sum())])
            if tr_t.sum() < 10:
                continue
            sc = StandardScaler().fit(Xtr)
            clf = LinearSVC(C=1.0, class_weight="balanced", dual="auto", max_iter=20000)
            clf.fit(sc.transform(Xtr), ytr)
            oof_t[fold_t == f] = clf.decision_function(sc.transform(Xt[fold_t == f]))
            oof_d[fold_d == f] = clf.decision_function(sc.transform(Xd[fold_d == f]))
        q = qvalues(oof_t, oof_d)
        new_pos = q <= a.q
        print(f"  iter {it}: positives {pos.sum()} -> {new_pos.sum()}  "
              f"(targets at q<= {a.q}: {new_pos.sum()})")
        if new_pos.sum() == pos.sum() and (new_pos == pos).all():
            pos = new_pos
            break
        pos = new_pos if new_pos.sum() >= 10 else pos

    # Diagnostic: which features the final model leans on (standardized coefficients). If the model
    # separates target-vs-decoy by envelope-shape features and near-ignores log_intensity, that is why
    # it underperforms an intensity ranking for real-vs-junk.
    Xtr = np.vstack([Xt[pos], Xd])
    ytr = np.concatenate([np.ones(pos.sum()), np.zeros(nd)])
    sc = StandardScaler().fit(Xtr)
    clf = LinearSVC(C=1.0, class_weight="balanced", dual="auto", max_iter=20000).fit(sc.transform(Xtr), ytr)
    tr_acc = (clf.predict(sc.transform(Xtr)) == ytr).mean()
    print(f"\nfinal model target-vs-decoy train accuracy: {100*tr_acc:.1f}%  (near 100% => fully separable)")
    print("SVM standardized weights (|w| desc):")
    for name, w in sorted(zip(FEATURE_NAMES, clf.coef_[0]), key=lambda x: -abs(x[1])):
        print(f"    {name:20s} {w:+.3f}")

    # ---- evaluation: SVM score vs intensity, recall-vs-retained-fraction ----
    if not a.gt:
        print("no --gt; skipping recall eval")
        return
    refs = load_refs(a.gt)
    n_ref = max(len(refs), 1)
    base_recall = recall_at(np.ones(nt, bool), meta_t, refs, a.rt)
    print(f"\nground-truth peaks: {len(refs)}   baseline recall (all {nt}): "
          f"{base_recall}/{len(refs)} = {100*base_recall/n_ref:.1f}%")
    print(f"\n{'frac':>5} {'kept':>8}  {'SVM recall':>10}  {'Inten recall':>12}  {'Decon recall':>12}")
    order_svm = np.argsort(-oof_t)
    order_int = np.argsort(-Xt[:, INTEN_IDX])
    order_dec = np.argsort(-Xt[:, DECON_IDX])
    for fr in [1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1]:
        k = max(1, int(round(nt * fr)))
        def rec(order):
            m = np.zeros(nt, bool); m[order[:k]] = True
            return 100 * recall_at(m, meta_t, refs, a.rt) / n_ref
        print(f"{fr:5.2f} {k:8d}  {rec(order_svm):9.1f}%  {rec(order_int):11.1f}%  {rec(order_dec):11.1f}%")


if __name__ == "__main__":
    main()
