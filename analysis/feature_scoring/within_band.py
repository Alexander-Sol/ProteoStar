#!/usr/bin/env python
"""Does envelope quality separate real from junk *within* an intensity band?

Tests the claim "low-intensity features are junk". Labels each target feature real/other by matching
to the PSM ground truth (mass +-20 ppm + RT window -- same as the recall scorer), then, within each
intensity quintile, asks whether shape scores (decon fit, ppm consistency, isotope count, charge
support, persistence) separate PSM-matched from unmatched features.

Two outputs:
  1. Match rate per intensity quintile -- how many real (PSM) features live in the low-intensity bins.
  2. Within-band AUC of each shape score for matched-vs-unmatched. AUC>0.5 => the score carries
     real-vs-other signal AT FIXED INTENSITY. (Conservative: "unmatched" is contaminated with real-
     but-unidentified features, so a clear AUC>0.5 is strong evidence; ~0.5 is ambiguous.)

Usage: within_band.py <target.tsv> <gt.tsv> <rt_delta>
"""
import csv, sys, bisect
import numpy as np

MASS_PPM = 20.0
SHAPE = ["decon", "neg_ppm", "max_iso", "num_cs", "persist"]


def load(path):
    rows = []
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        try:
            inten = float(r["Summed Intensity"])
            decon = float(r["Decon Score"])
            ppm = float(r["PPM Spread"]); ppm = 15.0 if ppm >= 999 else ppm
            iso = float(r["Max Num Isotopes"])
            ncs = float(r["Num Charge States"])
            persist = float(r["RT End"]) - float(r["RT Start"])
            mass = float(r["Monoisotopic Mass"]); rt = float(r["RT Apex"])
        except (KeyError, ValueError):
            continue
        rows.append(dict(inten=inten, decon=decon, neg_ppm=-ppm, max_iso=iso, num_cs=ncs,
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


def auc(pos, neg):
    """P(score_pos > score_neg) via Mann-Whitney; 0.5 = no separation."""
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    allv = np.concatenate([pos, neg])
    order = allv.argsort()
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(allv) + 1)
    # average ranks for ties
    _, inv, counts = np.unique(allv, return_inverse=True, return_counts=True)
    sums = np.zeros(len(counts)); np.add.at(sums, inv, ranks)
    ranks = (sums / counts)[inv]
    r_pos = ranks[: len(pos)].sum()
    return (r_pos - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def main():
    rows = load(sys.argv[1]); refs = load_refs(sys.argv[2]); rt_delta = float(sys.argv[3])
    n = len(rows)
    # label: feature matches a PSM peak (mass+RT)
    feats = sorted(rows, key=lambda x: x["mass"])
    masses = [f["mass"] for f in feats]
    for f in feats:
        f["real"] = False
    for (rmass, rrt) in refs:
        tol = rmass * MASS_PPM / 1e6
        lo = bisect.bisect_left(masses, rmass - tol); hi = bisect.bisect_right(masses, rmass + tol)
        for i in range(lo, hi):
            if abs(feats[i]["rt"] - rrt) <= rt_delta:
                feats[i]["real"] = True
                break
    inten = np.array([f["inten"] for f in feats])
    real = np.array([f["real"] for f in feats])
    print(f"features {n}   PSM-matched {int(real.sum())} ({100*real.mean():.2f}%)   refs {len(refs)}")

    # quintiles by intensity
    q = np.quantile(inten, [0.2, 0.4, 0.6, 0.8])
    band = np.digitize(inten, q)  # 0..4 low->high
    print(f"\nintensity quintile   n     matched   match-rate    frac-of-all-matched")
    tot_matched = int(real.sum())
    for b in range(5):
        m = band == b
        nm = int(real[m].sum())
        print(f"  Q{b+1} ({'low' if b==0 else 'high' if b==4 else '  '})        {m.sum():6d}  {nm:6d}    "
              f"{100*nm/max(m.sum(),1):6.2f}%       {100*nm/max(tot_matched,1):6.1f}%")

    print(f"\nwithin-band AUC (matched vs unmatched), per shape score:")
    print(f"  {'band':>6} {'n_matched':>9} " + " ".join(f"{s:>8}" for s in SHAPE))
    for b in range(5):
        m = band == b
        idx = np.nonzero(m)[0]
        rb = real[idx]
        line = f"  Q{b+1:<5} {int(rb.sum()):9d} "
        for s in SHAPE:
            vals = np.array([feats[i][s] for i in idx])
            line += f" {auc(vals[rb], vals[~rb]):8.3f}"
        print(line)
    print("\n(AUC>0.5 => that shape score ranks real above unmatched at fixed intensity; "
          "conservative since 'unmatched' contains real-but-unidentified features.)")


if __name__ == "__main__":
    main()
