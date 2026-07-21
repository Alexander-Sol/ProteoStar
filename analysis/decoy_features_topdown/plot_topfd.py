#!/usr/bin/env python3
"""Score distribution of TopFD's Jurkat features under OUR envelope-fit cosine, vs our own features.

TopFD (an established top-down feature detector) features were parsed from its *_feature.xml, expanded to
one row per (proteoform, charge), and scored by decoy_score_export against averagine 1.0 on the same
Jurkat raw and with the same settings we score our own features. This overlays TopFD's 2-D and 3-D score
distribution against our target (ID-matched + unmatched) and the two strong decoys, so we can see where an
independent tool's features land on the same axis.

Usage:  python plot_topfd.py
"""
import os, bisect, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "plots")
GT = os.path.join(HERE, "..", "topdown_bench", "gt_jurkat_intersect.tsv")
MASS_PPM = 20.0; C13 = 1.0033548; KMAX = 3; RT_TOL = 1.0
COLORS = {"topfd": "#7B3294", "matched": "#111111", "unmatched": "#888888",
          "weird": "#D55E00", "shuffled": "#009E73"}
SCORES = [("score2d", "2-D (apex) envelope-fit cosine"), ("score3d", "3-D (m/z x RT) envelope-fit cosine")]


def load(name):
    p = os.path.join(HERE, f"s_jurkat_{name}.tsv")
    return np.atleast_1d(np.genfromtxt(p, delimiter="\t", names=True))


def load_gt():
    g = np.atleast_1d(np.genfromtxt(GT, delimiter="\t", names=True))
    o = np.argsort(g["mono_mass"]); return g["mono_mass"][o], g["rt"][o]


def match_mask(feat, gm, gr):
    out = np.zeros(len(feat), bool)
    for i in range(len(feat)):
        m = feat["mono"][i]; rt = feat["rt"][i]
        for k in range(-KMAX, KMAX + 1):
            c = m + k * C13; t = c * MASS_PPM * 1e-6
            lo = bisect.bisect_left(gm, c - t); hi = bisect.bisect_right(gm, c + t)
            if lo < hi and np.any(np.abs(gr[lo:hi] - rt) <= RT_TOL):
                out[i] = True; break
    return out


def main():
    gm, gr = load_gt()
    topfd = load("topfd"); tgt = load("target")
    weird = load("weird"); shuf = load("shuffled")
    tm = match_mask(tgt, gm, gr)
    fm = match_mask(topfd, gm, gr)
    print(f"TopFD  n={len(topfd)}  ID-matched={fm.sum()} ({100*fm.mean():.1f}%)  "
          f"med 2D={np.median(topfd['score2d']):.3f}  med 3D={np.median(topfd['score3d']):.3f}")
    print(f"ours   target ID-matched med 2D={np.median(tgt['score2d'][tm]):.3f}  n_matched={tm.sum()}")

    # decoy-calibrated cut: Jurkat weird 95%-decoy-rejection 2-D threshold (5% of weird decoys pass)
    thr = None
    js = os.path.join(HERE, "topdown_decoy_summary.json")
    if os.path.exists(js):
        d = json.load(open(js))
        thr = d.get("jurkat", {}).get("decoys", {}).get("weird", {}).get("score2d", {}).get("prec95", {}).get("thr")
    if thr:
        print(f"\nat the weird-decoy 5%-pass 2-D threshold ({thr:.3f}):")
        print(f"  TopFD features passing:        {100*np.mean(topfd['score2d']>=thr):.1f}%")
        print(f"  our ID-matched targets passing:{100*np.mean(tgt['score2d'][tm]>=thr):.1f}%")

    fig, axes = plt.subplots(1, 2, figsize=(13, 6.4))
    bins = np.linspace(0, 1, 60)
    for ax, (score, xlabel) in zip(axes, SCORES):
        ax.hist(tgt[score][tm], bins=bins, density=True, histtype="step", lw=1.8,
                color=COLORS["matched"], label=f"ours: target, ID-matched (n={tm.sum()})")
        ax.hist(tgt[score][~tm], bins=bins, density=True, histtype="step", lw=1.2, ls=":",
                color=COLORS["unmatched"], label=f"ours: target, unmatched (n={(~tm).sum()})")
        ax.hist(weird[score], bins=bins, density=True, histtype="step", lw=1.4,
                color=COLORS["weird"], label=f"ours: weird decoy (n={len(weird)})")
        ax.hist(shuf[score], bins=bins, density=True, histtype="step", lw=1.4,
                color=COLORS["shuffled"], label=f"ours: shuffled decoy (n={len(shuf)})")
        ax.hist(topfd[score], bins=bins, density=True, histtype="stepfilled", lw=2.2, alpha=0.35,
                edgecolor=COLORS["topfd"], facecolor=COLORS["topfd"],
                label=f"TopFD, all (n={len(topfd)})")
        if thr and score == "score2d":
            ax.axvline(thr, color="#444444", lw=1.0, ls="--")
            ax.text(thr + 0.01, ax.get_ylim()[1]*0.9, "weird decoy\n5%-pass cut", fontsize=7, color="#444444")
        ax.set_xlabel(xlabel); ax.set_ylabel("density")
        ax.set_title(score.replace("score", "").upper() + " score")
        ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Jurkat — TopFD features vs ours, scored by the same envelope-fit cosine",
                 fontsize=13, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(OUT, "topfd_jurkat.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("\nwrote plots/topfd_jurkat.png")


if __name__ == "__main__":
    main()
