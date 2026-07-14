#!/usr/bin/env python3
"""ROC curves for target/decoy separation by envelope-fit cosine.

Positive class = ID-matched target features; negative class = each decoy set. One ROC per decoy, drawn
for the 2-D (apex) and 3-D (m/z x RT) score, per benchmark file. AUC is the area under each curve
(equals the Mann-Whitney statistic reported in plot_decoys.py). Also saves the FPR/TPR arrays and AUCs
to roc_data.json for the report.

Usage:  python roc_curves.py
"""
import os, bisect, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "plots")
os.makedirs(OUT, exist_ok=True)

FILES = [
    ("10", r"F:\flashlfq-rust\analysis\dinosaur_bench\gt_10min.tsv", 0.5, "10-min (CA/Lumos)"),
    ("65", r"F:\flashlfq-rust\analysis\longer_gradients\medium_ground_truth.tsv", 1.0, "65-min (glyco)"),
    ("2h", r"F:\flashlfq-rust\analysis\longer_gradients\long_ground_truth_msms.tsv", 0.5, "2-hr (IonStar)"),
]
DECOYS = ["shifted", "weird", "shuffled"]
LABELS = {"shifted": "shifted spacing (0.94 Da)", "weird": "weird averagine (Fe+Cl)", "shuffled": "shuffled envelope"}
COLORS = {"shifted": "#0072B2", "weird": "#D55E00", "shuffled": "#009E73"}
MASS_PPM = 20.0

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

def roc(pos, neg):
    """FPR, TPR sweeping threshold high->low; AUC by trapezoid."""
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    s = np.concatenate([pos, neg])
    order = np.argsort(-s, kind="mergesort")
    y = y[order]
    tp = np.cumsum(y); fp = np.cumsum(1 - y)
    tpr = np.concatenate([[0.0], tp / max(len(pos), 1)])
    fpr = np.concatenate([[0.0], fp / max(len(neg), 1)])
    auc = float(np.trapezoid(tpr, fpr))
    return fpr, tpr, auc

def main():
    data = {}
    for fk, gt_path, rt_tol, title in FILES:
        tgt = load(fk, "target"); gm, gr = load_gt(gt_path)
        matched = match_mask(tgt, gm, gr, rt_tol)
        data[fk] = {}
        fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
        for ax, score, sl in zip(axes, ("score2d", "score3d"), ("2-D (apex)", "3-D (m/z x RT)")):
            pos = tgt[score][matched]
            for d in DECOYS:
                neg = load(fk, d)[score]
                fpr, tpr, a = roc(pos, neg)
                # subsample points for a light figure / json
                if len(fpr) > 2000:
                    sel = np.linspace(0, len(fpr) - 1, 2000).astype(int)
                    fpr_s, tpr_s = fpr[sel], tpr[sel]
                else:
                    fpr_s, tpr_s = fpr, tpr
                ax.plot(fpr_s, tpr_s, color=COLORS[d], lw=2.0, label=f"{LABELS[d]}  AUC={a:.3f}")
                data[fk].setdefault(d, {})[score] = a
            ax.plot([0, 1], [0, 1], color="#999999", lw=1.0, ls=":")
            ax.set_xlabel("false-positive rate (decoys accepted)")
            ax.set_ylabel("true-positive rate (matched targets accepted)")
            ax.set_title(f"{sl} envelope-fit cosine")
            ax.set_xlim(0, 1); ax.set_ylim(0, 1.001)
            ax.legend(fontsize=8, loc="lower right")
            ax.spines[["top", "right"]].set_visible(False)
        fig.suptitle(f"{title} — ROC: ID-matched target vs decoy", fontsize=13, weight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        p = os.path.join(OUT, f"roc_{fk}.png")
        fig.savefig(p, dpi=140); plt.close(fig)
        print(f"[{fk}] wrote {p}")
    with open(os.path.join(HERE, "roc_data.json"), "w") as fh:
        json.dump(data, fh, indent=2)
    print("wrote roc_data.json")

if __name__ == "__main__":
    main()
