#!/usr/bin/env python
"""Test the isotopologue co-elution correlation feature (IsoCorr All/Top5/Top3).

Question: does isotope-trace correlation separate real (PSM-matched) from junk better than the
envelope-fit scores did -- especially on the large files where Decon/PPM carried ~no signal (AUC ~0.5)?

Reports, for each score column, the overall AUC (matched vs unmatched) and the per-intensity-quintile
AUC. AUC>0.5 => the score ranks real above unmatched. Decon Score is included as the fit-score
reference. IsoCorr sentinel (-2, uncomputable: <2 teeth or <3 scans) is kept as-is: it sorts as
"low correlation" = junk-like, which is the right treatment.

Usage: iso_corr_test.py <target.tsv> <gt.tsv> <rt_delta>
"""
import csv, sys, bisect
import numpy as np

MASS_PPM = 20.0
SCORES = ["IsoCorr All", "IsoCorr Top5", "IsoCorr Top3", "Decon Score", "PPM Spread"]


def load(path):
    rows = []
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        try:
            d = {"log_int": np.log1p(float(r["Summed Intensity"])),
                 "mass": float(r["Monoisotopic Mass"]), "rt": float(r["RT Apex"])}
            for s in SCORES:
                v = float(r[s])
                # PPM Spread: lower is better -> negate so "higher = better" holds for all columns.
                d[s] = (-v if v < 999 else -20.0) if s == "PPM Spread" else v
            rows.append(d)
        except (KeyError, ValueError):
            continue
    return rows


def load_refs(path):
    refs = []
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        try:
            refs.append((float(r["mono_mass"]), float(r["rt"])))
        except (ValueError, KeyError):
            continue
    return refs


def auc(pos, neg):
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    allv = np.concatenate([pos, neg])
    r = np.argsort(np.argsort(allv)) + 1.0
    # average ranks for ties
    _, inv, counts = np.unique(allv, return_inverse=True, return_counts=True)
    sums = np.zeros(len(counts)); np.add.at(sums, inv, r); r = (sums / counts)[inv]
    return (r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def main():
    rows = load(sys.argv[1]); refs = load_refs(sys.argv[2]); rt_delta = float(sys.argv[3])
    order = np.argsort([x["mass"] for x in rows]); masses = [rows[i]["mass"] for i in order]
    y = np.zeros(len(rows), bool)
    for (rmass, rrt) in refs:
        tol = rmass * MASS_PPM / 1e6
        lo = bisect.bisect_left(masses, rmass - tol); hi = bisect.bisect_right(masses, rmass + tol)
        for j in range(lo, hi):
            if abs(rows[order[j]]["rt"] - rrt) <= rt_delta:
                y[order[j]] = True
                break
    li = np.array([x["log_int"] for x in rows])
    edges = np.quantile(li, [.2, .4, .6, .8]); band = np.digitize(li, edges)
    print(f"features {len(rows)}  PSM-matched {int(y.sum())}  refs {len(refs)}")
    # sentinel fraction for the iso columns
    for s in ["IsoCorr All", "IsoCorr Top5", "IsoCorr Top3"]:
        v = np.array([x[s] for x in rows]); print(f"  {s}: {100*np.mean(v<=-2):.1f}% uncomputable, median(computable)={np.median(v[v>-2]):.3f}")
    print(f"\n{'score':>14} {'overall':>8}   " + " ".join(f"Q{b+1:>4}" for b in range(5)))
    for s in SCORES:
        v = np.array([x[s] for x in rows])
        row = f"{s:>14} {auc(v[y], v[~y]):8.3f}   "
        for b in range(5):
            m = band == b
            row += f"{auc(v[m & y], v[m & ~y]):5.2f} "
        print(row)
    print("\n(AUC>0.5 => ranks real above unmatched; per-quintile = at fixed intensity.)")


if __name__ == "__main__":
    main()
