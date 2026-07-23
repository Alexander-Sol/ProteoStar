#!/usr/bin/env python
"""Test the intensity x fit(+persistence) 2-D gate: does combining intensity with envelope-quality
features rank features better than intensity ALONE for retaining PSM recall?

Supervised logistic regression, label = PSM-matched (mass +-20 ppm + RT window). Evaluated OUT-OF-FOLD
(each feature scored by a model that never saw it) so it is a fair generalization test, not the
circular train==test. Compares nested feature sets so we can read what each term adds:
  int            : log_intensity only  (== the intensity ranking baseline)
  int+decon      : + envelope-fit cosine
  int+decon+pers : + chromatographic persistence (RT extent)
  full           : + min_decon, ppm_spread(+missing), max_iso, num_cs, cross_charge, num_members

CAVEAT baked into the interpretation: the PSM label is abundance-biased (an ID needs MS2 signal), so
this can only measure "retains identified peptides better", NOT rescue of genuine low-abundance
features. A win here is real but a floor on the true value; a null is not proof of no value.

Usage: combo_gate.py <target.tsv> <gt.tsv> <rt_delta>
"""
import csv, sys, bisect
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

MASS_PPM = 20.0
ALL_FEATS = ["log_intensity", "decon", "persist", "min_decon", "ppm", "ppm_missing",
             "max_iso", "num_cs", "cross_charge", "num_members"]
MODELS = {
    "int": ["log_intensity"],
    "int+isocorr": ["log_intensity", "isocorr"],
    "int+decon": ["log_intensity", "decon"],
    "full": ALL_FEATS,
    "full+isocorr": ALL_FEATS + ["isocorr"],
}


def load(path):
    rows = []
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        try:
            inten = float(r["Summed Intensity"]); decon = float(r["Decon Score"])
            mindec = float(r["Min Decon Score"]); ppm = float(r["PPM Spread"])
            iso = float(r["Max Num Isotopes"]); ncs = float(r["Num Charge States"])
            cc = float(r["Cross-Charge Support"]); nm = float(r["Num Members"])
            persist = float(r["RT End"]) - float(r["RT Start"])
            mass = float(r["Monoisotopic Mass"]); rt = float(r["RT Apex"])
        except (KeyError, ValueError):
            continue
        rows.append(dict(log_intensity=np.log1p(inten), decon=decon, min_decon=mindec,
                         isocorr=float(r.get("IsoCorr Top3", -2.0)),
                         ppm=15.0 if ppm >= 999 else ppm, ppm_missing=1.0 if ppm >= 999 else 0.0,
                         max_iso=iso, num_cs=ncs, cross_charge=cc, num_members=nm,
                         persist=persist, mass=mass, rt=rt))
    return rows


def load_refs(path):
    refs = []
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        try:
            refs.append((float(r["mono_mass"]), float(r["rt"])))
        except (ValueError, KeyError):
            continue
    return refs


def label(rows, refs, rt_delta):
    feats = sorted(rows, key=lambda x: x["mass"]); masses = [f["mass"] for f in feats]
    y = np.zeros(len(feats), bool)
    for (rmass, rrt) in refs:
        tol = rmass * MASS_PPM / 1e6
        lo = bisect.bisect_left(masses, rmass - tol); hi = bisect.bisect_right(masses, rmass + tol)
        for i in range(lo, hi):
            if abs(feats[i]["rt"] - rrt) <= rt_delta:
                y[i] = True
                break
    return feats, y


def recall_of_topk(order, k, feats, refs, rt_delta):
    keep = set(order[:k].tolist())
    sub = sorted((feats[i] for i in keep), key=lambda x: x["mass"])
    masses = [f["mass"] for f in sub]
    matched = 0
    for (rmass, rrt) in refs:
        tol = rmass * MASS_PPM / 1e6
        lo = bisect.bisect_left(masses, rmass - tol); hi = bisect.bisect_right(masses, rmass + tol)
        for i in range(lo, hi):
            if abs(sub[i]["rt"] - rrt) <= rt_delta:
                matched += 1
                break
    return matched


def oof_scores(feats, y, cols, folds=3):
    X = np.array([[f[c] for c in cols] for f in feats], float)
    rng = np.random.default_rng(0); fold = rng.integers(0, folds, len(feats))
    s = np.zeros(len(feats))
    coefs = []
    for k in range(folds):
        tr = fold != k
        sc = StandardScaler().fit(X[tr])
        clf = LogisticRegression(class_weight="balanced", max_iter=2000, C=1.0)
        clf.fit(sc.transform(X[tr]), y[tr])
        s[fold == k] = clf.decision_function(sc.transform(X[fold == k]))
        coefs.append(dict(zip(cols, clf.coef_[0])))
    mean_coef = {c: np.mean([cf[c] for cf in coefs]) for c in cols}
    return s, mean_coef


def main():
    rows = load(sys.argv[1]); refs = load_refs(sys.argv[2]); rt_delta = float(sys.argv[3])
    feats, y = label(rows, refs, rt_delta)
    n = len(feats); n_ref = max(len(refs), 1)
    base = recall_of_topk(np.arange(n), n, feats, refs, rt_delta)
    print(f"features {n}  PSM-matched {int(y.sum())}  refs {len(refs)}  baseline recall {100*base/n_ref:.1f}%")

    scores = {}
    for name, cols in MODELS.items():
        s, coef = oof_scores(feats, y, cols)
        scores[name] = s
        if name == "full":
            print("  full-model mean CV coefficients (standardized, desc):")
            for c, w in sorted(coef.items(), key=lambda x: -abs(x[1])):
                print(f"      {c:16s} {w:+.3f}")
    # raw intensity ranking is identical to model 'int' up to monotone transform; keep explicit.
    inten_raw = np.array([f["log_intensity"] for f in feats])

    fracs = [1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1]
    print(f"\n  recall(%) at retained top-fraction  (rt +-{rt_delta})")
    print(f"  {'frac':>5} {'kept':>8} " + " ".join(f"{m:>14}" for m in MODELS) + f" {'intensity_raw':>13}")
    ord_raw = np.argsort(-inten_raw)
    ords = {m: np.argsort(-scores[m]) for m in MODELS}
    for fr in fracs:
        k = max(1, int(round(n * fr)))
        cells = " ".join(f"{100*recall_of_topk(ords[m], k, feats, refs, rt_delta)/n_ref:13.1f}%" for m in MODELS)
        raw = 100 * recall_of_topk(ord_raw, k, feats, refs, rt_delta) / n_ref
        print(f"  {fr:5.2f} {k:8d} {cells} {raw:12.1f}%")


if __name__ == "__main__":
    main()
