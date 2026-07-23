#!/usr/bin/env python
"""Recall-vs-retained-fraction tradeoff for candidate feature scores.

The untargeted detector emits ~5-20x as many features as comparable tools (Dinosaur).
Goal: find a per-feature score that ranks the real features above the junk tail, so we
can keep the top X% and drop the rest without losing PSM recall.

For each requested score COLUMN, rank features by that score (descending = "keep high"),
retain the top fraction, and measure recall against the PSM ground truth using the SAME
mass +-20 ppm / RT window logic as score_any.py. Prints, per score, a curve of
(fraction retained, feature count, recall) so scores can be compared apples-to-apples:
the best score holds recall highest as the retained fraction shrinks.

A score column may be an existing resolved-TSV column (e.g. "Summed Intensity",
"Num Charge States", "Cross-Charge Support") or one added later. "Summed Intensity" is
the intensity/noise-floor baseline to beat.

Usage:
  score_tradeoff.py <resolved.tsv> <ground_truth.tsv> <rt_delta> <score_col>[,<score_col>...]
                    [--fracs=1,0.9,...] [--asc=col1,col2]  (score cols where LOW is better)
"""
import csv, sys, bisect

MASS_PPM = 20.0
DEFAULT_FRACS = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05]


def load_features(path, score_cols):
    """Return list of (mass, rt, set(charges), {col: float_value}). Rows missing any
    requested score column (unparseable) are kept with that score as None and excluded
    from that score's ranking."""
    feats = []
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                mass = float(row["Monoisotopic Mass"])
                rt = float(row["RT Apex"])
            except (ValueError, KeyError):
                continue
            charges = set()
            for c in row.get("Charge States", "").split(";"):
                c = c.strip()
                if c:
                    try:
                        charges.add(int(c))
                    except ValueError:
                        pass
            scores = {}
            for col in score_cols:
                try:
                    scores[col] = float(row[col])
                except (ValueError, KeyError):
                    scores[col] = None
            feats.append((mass, rt, charges, scores))
    return feats


def load_refs(path):
    refs = []
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                mass = float(row["mono_mass"]); rt = float(row["rt"]); z = int(row["charge"])
            except (ValueError, KeyError):
                continue
            refs.append((mass, rt, z))
    return refs


def recall(feats, refs, rt_delta):
    """feats: list of (mass, rt, charges). Returns (matched, matched_charge)."""
    fs = sorted(feats, key=lambda x: x[0])
    masses = [f[0] for f in fs]
    matched = matched_charge = 0
    for (rmass, rrt, rz) in refs:
        tol = rmass * MASS_PPM / 1e6
        lo = bisect.bisect_left(masses, rmass - tol)
        hi = bisect.bisect_right(masses, rmass + tol)
        for i in range(lo, hi):
            if abs(fs[i][1] - rrt) <= rt_delta:
                matched += 1
                if rz in fs[i][2]:
                    matched_charge += 1
                break
    return matched, matched_charge


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = {a.split("=", 1)[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    fpath, rpath, rt_delta = args[0], args[1], float(args[2])
    score_cols = args[3].split(",")
    fracs = [float(x) for x in opts.get("--fracs", "").split(",") if x] or DEFAULT_FRACS
    asc = set(opts.get("--asc", "").split(",")) - {""}

    feats = load_features(fpath, score_cols)
    refs = load_refs(rpath)
    n_ref = max(len(refs), 1)
    base_m, base_mc = recall([(f[0], f[1], f[2]) for f in feats], refs, rt_delta)
    print(f"features: {len(feats)}   ref peaks: {len(refs)}   rt_delta {rt_delta}")
    print(f"baseline (all features): recall {base_m}/{len(refs)} = {100*base_m/n_ref:.1f}%   "
          f"charge {100*base_mc/n_ref:.1f}%\n")

    for col in score_cols:
        ranked = [f for f in feats if f[3][col] is not None]
        # descending = keep-high; ascending (--asc) = keep-low (e.g. ppm spread).
        ranked.sort(key=lambda f: f[3][col], reverse=(col not in asc))
        direction = "low=better" if col in asc else "high=better"
        print(f"=== score: {col}  ({direction}, {len(ranked)} scored) ===")
        print(f"  {'frac':>5} {'kept':>8} {'recall':>8} {'charge':>7}  {'recall/base':>11}")
        for fr in fracs:
            k = max(1, int(round(len(ranked) * fr)))
            kept = ranked[:k]
            m, mc = recall([(f[0], f[1], f[2]) for f in kept], refs, rt_delta)
            print(f"  {fr:5.2f} {k:8d} {100*m/n_ref:7.1f}% {100*mc/n_ref:6.1f}% "
                  f"{100*m/max(base_m,1):10.1f}%")
        print()


if __name__ == "__main__":
    main()
