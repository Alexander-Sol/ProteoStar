#!/usr/bin/env python
"""Label-free fixed blend: rank by z(log_intensity) + lam*z(decon) -- no training, no decoy.

Decon Score is a detector output, so this needs no labels. Tests whether the modest per-band signal
that band_signal_check found (raw decon AUC ~0.6 on the large files) is captured by a simple fixed
blend, vs the failed learned within-band score. PSM GT used only for eval.

Usage: fixed_blend.py <target.tsv> <gt.tsv> <rt_delta> [lams]
"""
import csv, sys, bisect
import numpy as np

MASS_PPM = 20.0


def load(path):
    rows = []
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        try:
            rows.append({"li": np.log1p(float(r["Summed Intensity"])), "decon": float(r["Decon Score"]),
                         "mass": float(r["Monoisotopic Mass"]), "rt": float(r["RT Apex"])})
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


def recall_topk(score, k, rows, refs, rt, n_ref):
    keep = np.argpartition(-score, k - 1)[:k] if k < len(score) else np.arange(len(score))
    sub = sorted((rows[i] for i in keep), key=lambda x: x["mass"]); masses = [f["mass"] for f in sub]
    m = 0
    for (rm, rr) in refs:
        tol = rm * MASS_PPM / 1e6
        lo = bisect.bisect_left(masses, rm - tol); hi = bisect.bisect_right(masses, rm + tol)
        for i in range(lo, hi):
            if abs(sub[i]["rt"] - rr) <= rt:
                m += 1
                break
    return 100 * m / n_ref


def main():
    rows = load(sys.argv[1]); refs = load_refs(sys.argv[2]); rt = float(sys.argv[3])
    lams = [float(x) for x in (sys.argv[4].split(",") if len(sys.argv) > 4 else "0,0.25,0.5,1".split(","))]
    n = len(rows); n_ref = max(len(refs), 1)
    li = np.array([r["li"] for r in rows]); dec = np.array([r["decon"] for r in rows])
    zi = (li - li.mean()) / li.std(); zd = (dec - dec.mean()) / dec.std()
    print(f"features {n}  refs {len(refs)}")
    print(f"  {'frac':>5} " + " ".join(f"lam={l:<4}" for l in lams) + f"  {'intensity':>9}")
    for fr in [0.4, 0.3, 0.25, 0.2, 0.15, 0.1]:
        k = max(1, int(round(n * fr)))
        cells = " ".join(f"{recall_topk(zi + l * zd, k, rows, refs, rt, n_ref):8.1f}" for l in lams)
        ri = recall_topk(zi, k, rows, refs, rt, n_ref)
        print(f"  {fr:5.2f} {cells}  {ri:9.1f}")


if __name__ == "__main__":
    main()
