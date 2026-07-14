#!/usr/bin/env python3
"""Decoy-vs-target score analysis for the untargeted MS1 feature detector.

For each benchmark file we detected features four ways:
  target    — the real averagine envelope (the normal pipeline)
  shifted   — 0.94-Da shifted isotope spacing (off-lattice decoy)
  weird     — exotic Fe+Cl averagine ("P2 C1 N1 Cl1 Fe1", composition replaces backbone)
  shuffled  — real averagine weights randomly permuted (on-lattice ratio decoy)

Each feature carries a 2-D (apex-scan m/z) and 3-D (m/z x RT window) envelope-fit cosine, scored
against the SAME envelope it was detected with. Target features are split into ID-matched (matching a
reference PSM/quantified peak in mass+RT) vs unmatched. We overlay the ID-matched-target score
distribution against each decoy's, per file, for both scores, and report ROC-AUC separation
(matched-target = positive, decoy = negative; AUC 0.5 = no separation, 1.0 = perfect).

Usage:  python plot_decoys.py
"""
import os
import bisect
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

# (file key, ground-truth table, RT match tolerance in minutes)
FILES = [
    ("10", r"F:\flashlfq-rust\analysis\dinosaur_bench\gt_10min.tsv", 0.5),
    ("65", r"F:\flashlfq-rust\analysis\longer_gradients\medium_ground_truth.tsv", 1.0),
    ("2h", r"F:\flashlfq-rust\analysis\longer_gradients\long_ground_truth_msms.tsv", 0.5),
]
FILE_TITLES = {"10": "10-min (CA/Lumos)", "65": "65-min (glyco)", "2h": "2-hr (IonStar)"}

DECOYS = ["shifted", "weird", "shuffled"]
DECOY_LABELS = {
    "shifted": "shifted spacing (0.94 Da)",
    "weird": "weird averagine (Fe+Cl)",
    "shuffled": "shuffled envelope",
}
# Okabe-Ito colourblind-safe palette
COLORS = {
    "matched":  "#111111",     # ID-matched target
    "unmatched": "#111111",   # unmatched target (context)
    "shifted": "#0072B2",     # blue
    "weird": "#D55E00",       # vermillion
    "shuffled": "#009E73",    # green
}
MASS_PPM = 20.0


def load_scores(fk, model):
    path = os.path.join(HERE, f"s{fk}_{model}.tsv")
    if not os.path.exists(path):
        return None
    rows = np.genfromtxt(path, delimiter="\t", names=True)
    if rows.size == 0:
        return None
    return np.atleast_1d(rows)


def load_gt(path):
    gt = np.genfromtxt(path, delimiter="\t", names=True)
    gt = np.atleast_1d(gt)
    order = np.argsort(gt["mono_mass"])
    return gt["mono_mass"][order], gt["rt"][order]


def match_mask(feat, gt_mass, gt_rt, rt_tol):
    """Boolean mask: which features match any GT peak within MASS_PPM (mass) and rt_tol (RT)."""
    masses = feat["mono"]
    rts = feat["rt"]
    out = np.zeros(len(masses), dtype=bool)
    for i in range(len(masses)):
        m = masses[i]
        tol = m * MASS_PPM * 1e-6
        lo = bisect.bisect_left(gt_mass, m - tol)
        hi = bisect.bisect_right(gt_mass, m + tol)
        if lo == hi:
            continue
        # any GT peak in the mass window also within RT tolerance?
        if np.any(np.abs(gt_rt[lo:hi] - rts[i]) <= rt_tol):
            out[i] = True
    return out


def auc(pos, neg):
    """ROC-AUC via Mann-Whitney U (pos = matched target scores, neg = decoy scores)."""
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    allv = np.concatenate([pos, neg])
    ranks = allv.argsort().argsort().astype(float) + 1.0
    # average ranks for ties
    order = np.argsort(allv)
    sv = allv[order]
    r = ranks[order]
    i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        if j > i:
            r[i:j + 1] = (r[i] + r[j]) / 2.0
        i = j + 1
    ranks[order] = r
    rank_pos = ranks[:len(pos)]
    u = rank_pos.sum() - len(pos) * (len(pos) + 1) / 2.0
    return u / (len(pos) * len(neg))


def main():
    outdir = os.path.join(HERE, "plots")
    os.makedirs(outdir, exist_ok=True)
    summary = []  # (file, decoy, score, auc)

    for fk, gt_path, rt_tol in FILES:
        tgt = load_scores(fk, "target")
        if tgt is None:
            print(f"[{fk}] target scores missing — skip")
            continue
        gt_mass, gt_rt = load_gt(gt_path)
        matched = match_mask(tgt, gt_mass, gt_rt, rt_tol)
        n_match = int(matched.sum())
        print(f"[{fk}] target n={len(tgt)}  ID-matched={n_match} ({100*n_match/len(tgt):.1f}%)")

        decoy_scores = {d: load_scores(fk, d) for d in DECOYS}

        fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
        for ax, score in zip(axes, ("score2d", "score3d")):
            tgt_m = tgt[score][matched]
            tgt_u = tgt[score][~matched]
            bins = np.linspace(0, 1, 60)
            # unmatched target (context) and matched target
            ax.hist(tgt_u, bins=bins, density=True, histtype="step", lw=2.0,
                    color=COLORS["unmatched"], label=f"target, unmatched (n={len(tgt_u)})")
            ax.hist(tgt_m, bins=bins, density=True, histtype="stepfilled", lw=2.0, alpha=0.35,
                    edgecolor=COLORS["matched"], facecolor=COLORS["matched"],
                    label=f"target, ID-matched (n={len(tgt_m)})")
            for d in DECOYS:
                ds = decoy_scores[d]
                if ds is None:
                    continue
                a = auc(tgt_m, ds[score])
                summary.append((fk, d, score, a))
                ax.hist(ds[score], bins=bins, density=True, histtype="step", lw=2.0,
                        color=COLORS[d], label=f"{DECOY_LABELS[d]} (n={len(ds)}, AUC={a:.3f})")
            ax.set_xlabel(f"{'2-D (apex)' if score=='score2d' else '3-D (m/z x RT)'} envelope-fit cosine")
            ax.set_ylabel("density")
            ax.set_title(f"{'2-D' if score=='score2d' else '3-D'} score")
            ax.legend(fontsize=8, loc="upper right")
            ax.spines[["top", "right"]].set_visible(False)
        fig.suptitle(f"{FILE_TITLES[fk]} — ID-matched target vs decoys", fontsize=13, weight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        outpng = os.path.join(outdir, f"decoy_scores_{fk}.png")
        fig.savefig(outpng, dpi=140)
        plt.close(fig)
        print(f"[{fk}] wrote {outpng}")

    # AUC summary table
    print("\n=== ROC-AUC (matched-target vs decoy; higher = better separation) ===")
    print(f"{'file':<5}{'decoy':<12}{'2-D':>8}{'3-D':>8}")
    sumpath = os.path.join(outdir, "auc_summary.tsv")
    with open(sumpath, "w") as fh:
        fh.write("file\tdecoy\tauc_2d\tauc_3d\n")
        by = {}
        for fk, d, score, a in summary:
            by.setdefault((fk, d), {})[score] = a
        for (fk, d), sc in by.items():
            a2, a3 = sc.get("score2d", float("nan")), sc.get("score3d", float("nan"))
            print(f"{fk:<5}{d:<12}{a2:>8.3f}{a3:>8.3f}")
            fh.write(f"{fk}\t{d}\t{a2:.4f}\t{a3:.4f}\n")
    print(f"\nwrote {sumpath}")


if __name__ == "__main__":
    main()
