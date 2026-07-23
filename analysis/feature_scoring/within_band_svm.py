#!/usr/bin/env python
"""Label-free within-band semi-supervised rescoring (option A).

Replaces combo_gate.py's PSM labels with a self-supervised scheme:
  - NEGATIVES: noise-decoy features (NOISE_DECOY_SHIFT run) -- real averagine comb shifted off the seed,
    so they are genuine noise coincidences (worse fit at the SAME intensity as targets; verified).
  - POSITIVES: seeded WITHIN each intensity band by apex-scan cosine (Decon Score) -- the best-fitting
    targets at each intensity level. Re-selected each iteration by the model's own out-of-fold score
    (semi-supervised). Doing it per-band keeps intensity from being the discriminator, forcing the model
    to learn a SHAPE-quality axis valid at every intensity.

The learned shape score carries no intensity signal by construction, so the final ranking BLENDS it with
intensity explicitly:  composite = z(log_intensity) + lam * z(shape_score). We sweep lam and compare the
recall-vs-retained-fraction curve to intensity alone (lam=0). PSM ground truth is used ONLY for
evaluation, never for training -- so a win here is label-free.

Usage: within_band_svm.py --target T.tsv --decoy D.tsv --gt GT.tsv --rt 0.5
                          [--bands 10] [--seed-q 0.5] [--iters 5] [--lams 0,0.5,1,2]
"""
import argparse, csv, bisect
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

MASS_PPM = 20.0
SHAPE = ["decon", "min_decon", "ppm", "ppm_missing", "max_iso", "num_cs", "cross_charge", "num_members", "persist"]


def load(path):
    rows = []
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        try:
            inten = float(r["Summed Intensity"])
            d = dict(log_int=np.log1p(inten), decon=float(r["Decon Score"]), min_decon=float(r["Min Decon Score"]),
                     max_iso=float(r["Max Num Isotopes"]), num_cs=float(r["Num Charge States"]),
                     cross_charge=float(r["Cross-Charge Support"]), num_members=float(r["Num Members"]),
                     persist=float(r["RT End"]) - float(r["RT Start"]),
                     mass=float(r["Monoisotopic Mass"]), rt=float(r["RT Apex"]))
            p = float(r["PPM Spread"]); d["ppm"] = 15.0 if p >= 999 else p; d["ppm_missing"] = 1.0 if p >= 999 else 0.0
        except (KeyError, ValueError):
            continue
        rows.append(d)
    return rows


def load_refs(path):
    refs = []
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        try:
            refs.append((float(r["mono_mass"]), float(r["rt"])))
        except (ValueError, KeyError):
            continue
    return refs


def psm_label(rows, refs, rt_delta):
    order = np.argsort([x["mass"] for x in rows])
    masses = [rows[i]["mass"] for i in order]
    y = np.zeros(len(rows), bool)
    for (rmass, rrt) in refs:
        tol = rmass * MASS_PPM / 1e6
        lo = bisect.bisect_left(masses, rmass - tol); hi = bisect.bisect_right(masses, rmass + tol)
        for j in range(lo, hi):
            i = order[j]
            if abs(rows[i]["rt"] - rrt) <= rt_delta:
                y[i] = True
                break
    return y


def recall_topk(score, k, rows, refs, rt_delta, n_ref):
    keep = np.argpartition(-score, k - 1)[:k] if k < len(score) else np.arange(len(score))
    sub = sorted((rows[i] for i in keep), key=lambda x: x["mass"])
    masses = [f["mass"] for f in sub]
    m = 0
    for (rmass, rrt) in refs:
        tol = rmass * MASS_PPM / 1e6
        lo = bisect.bisect_left(masses, rmass - tol); hi = bisect.bisect_right(masses, rmass + tol)
        for i in range(lo, hi):
            if abs(sub[i]["rt"] - rrt) <= rt_delta:
                m += 1
                break
    return 100 * m / n_ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True); ap.add_argument("--decoy", required=True)
    ap.add_argument("--gt", required=True); ap.add_argument("--rt", type=float, default=0.5)
    ap.add_argument("--bands", type=int, default=10); ap.add_argument("--seed-q", type=float, default=0.5)
    ap.add_argument("--iters", type=int, default=5); ap.add_argument("--lams", default="0,0.5,1,2,4")
    a = ap.parse_args()

    T = load(a.target); D = load(a.decoy); refs = load_refs(a.gt); n_ref = max(len(refs), 1)
    nt, nd = len(T), len(D)
    Xt_shape = np.array([[t[c] for c in SHAPE] for t in T]); Xd_shape = np.array([[d[c] for c in SHAPE] for d in D])
    li_t = np.array([t["log_int"] for t in T])
    # intensity bands from target quantiles; assign decoys by the same edges
    edges = np.quantile(li_t, np.linspace(0, 1, a.bands + 1)[1:-1])
    band_t = np.digitize(li_t, edges); band_d = np.digitize(np.array([d["log_int"] for d in D]), edges)
    ypsm = psm_label(T, refs, a.rt)
    print(f"targets {nt}  decoys {nd}  bands {a.bands}  PSM-matched {int(ypsm.sum())}  refs {len(refs)}")

    # seed positives: per band, top seed-q of targets by decon
    decon = Xt_shape[:, SHAPE.index("decon")]
    pos = np.zeros(nt, bool)
    for b in range(a.bands):
        m = band_t == b
        if m.sum() == 0:
            continue
        thr = np.quantile(decon[m], 1 - a.seed_q)
        pos |= m & (decon >= thr)

    rng = np.random.default_rng(0); fold = rng.integers(0, 3, nt)
    oof = np.zeros(nt)
    for it in range(a.iters):
        oof[:] = 0
        for f in range(3):
            tr = (fold != f) & pos
            if tr.sum() < 20:
                continue
            Xtr = np.vstack([Xt_shape[tr], Xd_shape]); ytr = np.concatenate([np.ones(tr.sum()), np.zeros(nd)])
            sc = StandardScaler().fit(Xtr)
            clf = LogisticRegression(class_weight="balanced", max_iter=2000).fit(sc.transform(Xtr), ytr)
            oof[fold == f] = clf.decision_function(sc.transform(Xt_shape[fold == f]))
        # re-select positives per band by oof score
        newpos = np.zeros(nt, bool)
        for b in range(a.bands):
            m = band_t == b
            if m.sum() == 0:
                continue
            thr = np.quantile(oof[m], 1 - a.seed_q)
            newpos |= m & (oof >= thr)
        changed = int((newpos != pos).sum()); pos = newpos
        print(f"  iter {it}: positives {int(pos.sum())}  (changed {changed})")

    # within-band AUC of the label-free shape score vs PSM label (eval only)
    def auc(s, y):
        pos_, neg_ = s[y], s[~y]
        if len(pos_) == 0 or len(neg_) == 0:
            return float("nan")
        r = np.argsort(np.argsort(np.concatenate([pos_, neg_]))) + 1
        return (r[:len(pos_)].sum() - len(pos_) * (len(pos_) + 1) / 2) / (len(pos_) * len(neg_))
    print("  shape-score within-band AUC vs PSM:", " ".join(
        f"Q{b+1}={auc(oof[band_t==b], ypsm[band_t==b]):.2f}" for b in range(a.bands)))

    # blend and evaluate: composite = z(log_int) + lam * z(shape)
    zi = (li_t - li_t.mean()) / li_t.std(); zs = (oof - oof.mean()) / oof.std()
    lams = [float(x) for x in a.lams.split(",")]
    base = recall_topk(np.ones(nt), nt, T, refs, a.rt, n_ref)
    print(f"\nbaseline recall (all): {base:.1f}%")
    fracs = [0.6, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1]
    print(f"  {'frac':>5} " + " ".join(f"lam={l:<4}" for l in lams) + f"  {'intens':>7} {'shape':>7}")
    for fr in fracs:
        k = max(1, int(round(nt * fr)))
        cells = " ".join(f"{recall_topk(zi + l * zs, k, T, refs, a.rt, n_ref):7.1f}" for l in lams)
        ri = recall_topk(zi, k, T, refs, a.rt, n_ref); rs = recall_topk(zs, k, T, refs, a.rt, n_ref)
        print(f"  {fr:5.2f} {cells}  {ri:7.1f} {rs:7.1f}")


if __name__ == "__main__":
    main()
