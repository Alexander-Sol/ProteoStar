#!/usr/bin/env python3
"""2-D cosine threshold at 95% recall of ID-matched targets, and how many non-ID (unmatched) target
features — and weird-averagine decoy features — clear that same threshold. Reuses the matching logic
from plot_decoys.py."""
import os, bisect
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FILES = [
    ("10", r"F:\flashlfq-rust\analysis\dinosaur_bench\gt_10min.tsv", 0.5, "10-min (CA/Lumos)"),
    ("65", r"F:\flashlfq-rust\analysis\longer_gradients\medium_ground_truth.tsv", 1.0, "65-min (glyco)"),
    ("2h", r"F:\flashlfq-rust\analysis\longer_gradients\long_ground_truth_msms.tsv", 0.5, "2-hr (IonStar)"),
]
MASS_PPM = 20.0
RECALL = 0.95

def load(fk, model):
    p = os.path.join(HERE, f"s{fk}_{model}.tsv")
    return np.atleast_1d(np.genfromtxt(p, delimiter="\t", names=True))

def load_gt(path):
    gt = np.atleast_1d(np.genfromtxt(path, delimiter="\t", names=True))
    o = np.argsort(gt["mono_mass"])
    return gt["mono_mass"][o], gt["rt"][o]

def match_mask(feat, gm, gr, rt_tol):
    out = np.zeros(len(feat), dtype=bool)
    for i in range(len(feat)):
        m = feat["mono"][i]; tol = m * MASS_PPM * 1e-6
        lo = bisect.bisect_left(gm, m - tol); hi = bisect.bisect_right(gm, m + tol)
        if lo < hi and np.any(np.abs(gr[lo:hi] - feat["rt"][i]) <= rt_tol):
            out[i] = True
    return out

print(f"{'file':<16}{'thr@95%rec':>12}{'matched':>10}{'unmatched>=thr':>16}{'unmatched tot':>15}"
      f"{'weird>=thr':>12}{'weird tot':>11}{'weird passrate':>16}")
for fk, gt_path, rt_tol, title in FILES:
    tgt = load(fk, "target"); weird = load(fk, "weird")
    gm, gr = load_gt(gt_path)
    matched = match_mask(tgt, gm, gr, rt_tol)
    s2 = tgt["score2d"]
    thr = float(np.percentile(s2[matched], (1 - RECALL) * 100))  # 5th pct of matched scores
    n_match = int(matched.sum())
    unmatched_scores = s2[~matched]
    n_unm_above = int(np.sum(unmatched_scores >= thr))
    n_unm_tot = int((~matched).sum())
    w2 = weird["score2d"]
    n_w_above = int(np.sum(w2 >= thr)); n_w_tot = len(w2)
    print(f"{title:<16}{thr:>12.4f}{n_match:>10d}{n_unm_above:>16d}{n_unm_tot:>15d}"
          f"{n_w_above:>12d}{n_w_tot:>11d}{100*n_w_above/n_w_tot:>15.2f}%")
