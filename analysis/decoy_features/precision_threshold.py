#!/usr/bin/env python3
"""2-D cosine threshold that REJECTS 95% of weird-averagine decoys (95th percentile of the decoy score
distribution — only 5% of decoys pass). Reports the resulting recall of ID-matched real features and how
many unmatched target features survive. Mirrors recall_threshold.py."""
import os, bisect
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FILES = [
    ("10", r"F:\flashlfq-rust\analysis\dinosaur_bench\gt_10min.tsv", 0.5, "10-min (CA/Lumos)"),
    ("65", r"F:\flashlfq-rust\analysis\longer_gradients\medium_ground_truth.tsv", 1.0, "65-min (glyco)"),
    ("2h", r"F:\flashlfq-rust\analysis\longer_gradients\long_ground_truth_msms.tsv", 0.5, "2-hr (IonStar)"),
]
MASS_PPM = 20.0
DECOY_REJECT = 0.95  # reject 95% of decoys -> threshold at 95th pct of decoy scores

def load(fk, model):
    return np.atleast_1d(np.genfromtxt(os.path.join(HERE, f"s{fk}_{model}.tsv"), delimiter="\t", names=True))

def load_gt(path):
    gt = np.atleast_1d(np.genfromtxt(path, delimiter="\t", names=True))
    o = np.argsort(gt["mono_mass"]); return gt["mono_mass"][o], gt["rt"][o]

def match_mask(feat, gm, gr, rt_tol):
    out = np.zeros(len(feat), dtype=bool)
    for i in range(len(feat)):
        m = feat["mono"][i]; tol = m * MASS_PPM * 1e-6
        lo = bisect.bisect_left(gm, m - tol); hi = bisect.bisect_right(gm, m + tol)
        if lo < hi and np.any(np.abs(gr[lo:hi] - feat["rt"][i]) <= rt_tol):
            out[i] = True
    return out

print(f"{'file':<16}{'thr@95%rej':>12}{'recall':>9}{'matched>=thr':>14}{'matched tot':>13}"
      f"{'unmatched>=thr':>16}{'unmatched tot':>15}{'decoys pass':>13}")
for fk, gt_path, rt_tol, title in FILES:
    tgt = load(fk, "target"); weird = load(fk, "weird")
    gm, gr = load_gt(gt_path)
    matched = match_mask(tgt, gm, gr, rt_tol)
    thr = float(np.percentile(weird["score2d"], DECOY_REJECT * 100))
    s2 = tgt["score2d"]
    m_scores = s2[matched]; u_scores = s2[~matched]
    n_m_above = int(np.sum(m_scores >= thr)); n_m_tot = int(matched.sum())
    n_u_above = int(np.sum(u_scores >= thr)); n_u_tot = int((~matched).sum())
    recall = n_m_above / n_m_tot
    dpass = float(np.mean(weird["score2d"] >= thr))
    print(f"{title:<16}{thr:>12.4f}{100*recall:>8.1f}%{n_m_above:>14d}{n_m_tot:>13d}"
          f"{n_u_above:>16d}{n_u_tot:>15d}{100*dpass:>12.2f}%")
