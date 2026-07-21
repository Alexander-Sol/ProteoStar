#!/usr/bin/env python
"""Top-down recall scorer with an off-by-one-aware diagnostic.

Same matching as analysis/dinosaur_bench/score_any.py (mass ppm + RT), but adds a second
"isotope-tolerant" pass that also counts a match when the detected mass is off from the
reference by k*(C13-C12) for k in [-KMAX, KMAX] within the ppm window. On large proteoforms
the monoisotope is many peaks below the envelope apex, so a ±1/±2 13C mono error is the
dominant, *fixable* failure mode -- separating it from a true miss tells us whether to spend
effort on mono assignment vs detection sensitivity.

Reports, per RT window:
    recall(strict)     -- mass within MASS_PPM, RT within delta        (== score_any)
    recall(isotope)    -- also allow k*1.00335 Da mono offset          (found the species)
    off-by-one gap     -- isotope - strict  (recall recoverable by fixing mono assignment)
    charge recall      -- ref charge in the matched feature's charge list (strict match)

Usage: score_topdown.py <features.tsv> <ground_truth.tsv> [rt_delta ...]
"""
import csv, sys, bisect

MASS_PPM = 20.0
C13 = 1.0033548
KMAX = 3


def load_features(path):
    feats = []
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        cols = set(r.fieldnames or [])
        is_dino = "rtApex" in cols and "charge" in cols
        for row in r:
            try:
                if is_dino:
                    mass = float(row["mass"]); rt = float(row["rtApex"])
                    charges = {int(float(row["charge"]))}
                else:
                    mass = float(row["Monoisotopic Mass"]); rt = float(row["RT Apex"])
                    charges = {int(c) for c in
                               (row.get("Charge States", "") or "").split(";") if c.strip()}
            except (ValueError, KeyError):
                continue
            feats.append((mass, rt, charges))
    feats.sort(key=lambda x: x[0])
    return feats


def load_refs(path):
    refs = []
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                refs.append((float(row["mono_mass"]), float(row["rt"]), int(float(row["charge"]))))
            except (ValueError, KeyError):
                continue
    return refs


def score(feats, refs, rt_delta):
    masses = [f[0] for f in feats]
    strict = strict_charge = iso = 0
    for (rmass, rrt, rz) in refs:
        # strict: single ppm window at k=0
        tol = rmass * MASS_PPM / 1e6
        lo = bisect.bisect_left(masses, rmass - tol)
        hi = bisect.bisect_right(masses, rmass + tol)
        hit = None
        for i in range(lo, hi):
            if abs(feats[i][1] - rrt) <= rt_delta:
                hit = feats[i]; break
        if hit is not None:
            strict += 1
            if rz in hit[2]:
                strict_charge += 1
            iso += 1
            continue
        # isotope-tolerant: try k*C13 offsets
        found = False
        for k in range(-KMAX, KMAX + 1):
            if k == 0:
                continue
            center = rmass + k * C13
            tol = center * MASS_PPM / 1e6
            lo = bisect.bisect_left(masses, center - tol)
            hi = bisect.bisect_right(masses, center + tol)
            for i in range(lo, hi):
                if abs(feats[i][1] - rrt) <= rt_delta:
                    found = True; break
            if found:
                break
        if found:
            iso += 1
    return strict, strict_charge, iso


def main():
    fpath, rpath = sys.argv[1], sys.argv[2]
    deltas = [float(x) for x in sys.argv[3:]] or [0.5]
    feats = load_features(fpath)
    refs = load_refs(rpath)
    n = max(len(refs), 1)
    print(f"features: {len(feats)}    ref peaks: {len(refs)}")
    for d in deltas:
        s, sc, iso = score(feats, refs, d)
        print(f"  RT +-{d:.2f} min:  strict {s}/{len(refs)}={100*s/n:.1f}%   "
              f"isotope-tol {iso}/{len(refs)}={100*iso/n:.1f}%   "
              f"off-by-one gap {100*(iso-s)/n:.1f}pp   "
              f"charge {sc}/{len(refs)}={100*sc/n:.1f}%")


if __name__ == "__main__":
    main()
