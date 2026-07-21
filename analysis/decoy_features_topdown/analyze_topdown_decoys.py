#!/usr/bin/env python3
"""Top-down decoy-vs-target score analysis.

For each top-down dataset we detected features four ways (envelope model is the only difference):
  target    - real averagine envelope
  shifted   - 0.94-Da shifted isotope spacing (off-lattice decoy)
  weird     - exotic Fe+Cl averagine (composition replaces the backbone)
  shuffled  - real averagine weights randomly permuted (on-lattice ratio decoy)

Every per-charge feature carries a 2-D (apex-scan m/z) and 3-D (m/z x RT window) envelope-fit cosine,
scored against the SAME envelope it was detected with. Target features are split into ID-matched
(matching an identified proteoform in mass+RT, isotope-tolerant) vs unmatched. We produce, per dataset:
  * score distributions: unmatched target, ID-matched target, and each decoy (2-D and 3-D);
  * ROC curves + AUC (positive = ID-matched target, negative = each decoy);
  * 95%-recall and 95%-precision operating thresholds per decoy strategy.

Isotope-tolerant matching mirrors score_topdown.py: a feature matches a proteoform if its mono is
within 20 ppm of the reference mono +/- k*(C13-C12) for k in [-3, 3], and its apex RT within RT_TOL.
This labels a feature as "detected the identified proteoform" even when the monoisotope is off-by-one
(the dominant top-down error), which is what we want for the positive class -- the envelope-fit score,
not mono precision, is the axis under test.

Usage:  python analyze_topdown_decoys.py
"""
import os, bisect, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "plots")
os.makedirs(OUT, exist_ok=True)

# (dataset key, ground-truth table, RT tolerance min, human title, n high-conf label)
DATASETS = [
    ("golden", os.path.join(HERE, "..", "topdown_bench", "gt_golden_all.tsv"), 1.0, "Golden (930 proteoforms)"),
    ("jurkat", os.path.join(HERE, "..", "topdown_bench", "gt_jurkat_intersect.tsv"), 1.0, "Jurkat (1541 hi-conf proteoforms)"),
]
DECOYS = ["shifted", "weird", "shuffled"]
LABELS = {"shifted": "shifted spacing (0.94 Da)", "weird": "weird averagine (Fe+Cl)", "shuffled": "shuffled envelope"}
COLORS = {"matched": "#111111", "unmatched": "#888888",
          "shifted": "#0072B2", "weird": "#D55E00", "shuffled": "#009E73"}
MASS_PPM = 20.0
C13 = 1.0033548
KMAX = 3
SCORES = [("score2d", "2-D (apex) envelope-fit cosine"), ("score3d", "3-D (m/z x RT) envelope-fit cosine")]


def load_scores(ds, model):
    path = os.path.join(HERE, f"s_{ds}_{model}.tsv")
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    rows = np.atleast_1d(np.genfromtxt(path, delimiter="\t", names=True))
    return rows if rows.size else None


def load_gt(path):
    gt = np.atleast_1d(np.genfromtxt(path, delimiter="\t", names=True))
    o = np.argsort(gt["mono_mass"])
    return gt["mono_mass"][o], gt["rt"][o]


def match_mask(feat, gm, gr, rt_tol):
    """Isotope-tolerant: feature matches a proteoform within 20 ppm at any k*C13 offset (k in +/-3) and RT."""
    out = np.zeros(len(feat), dtype=bool)
    masses = feat["mono"]; rts = feat["rt"]
    for i in range(len(masses)):
        m = masses[i]; rt = rts[i]
        for k in range(-KMAX, KMAX + 1):
            c = m + k * C13
            tol = c * MASS_PPM * 1e-6
            lo = bisect.bisect_left(gm, c - tol); hi = bisect.bisect_right(gm, c + tol)
            if lo < hi and np.any(np.abs(gr[lo:hi] - rt) <= rt_tol):
                out[i] = True
                break
    return out


def roc(pos, neg):
    """FPR, TPR (threshold high->low) and AUC by trapezoid. pos/neg are score arrays."""
    if len(pos) == 0 or len(neg) == 0:
        return np.array([0, 1.0]), np.array([0, 1.0]), float("nan")
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    s = np.concatenate([pos, neg])
    order = np.argsort(-s, kind="mergesort")
    y = y[order]
    tp = np.cumsum(y); fp = np.cumsum(1 - y)
    tpr = np.concatenate([[0.0], tp / len(pos)])
    fpr = np.concatenate([[0.0], fp / len(neg)])
    return fpr, tpr, float(np.trapezoid(tpr, fpr))


def op_points(pos, neg):
    """Two operating points, matching the bottom-up decoy convention (precision_threshold.py /
    recall_threshold.py):

    - recall95: threshold = 5th percentile of ID-matched target scores, so 95% of real IDs are
      retained. Report the decoy pass-rate (fraction of this decoy's features clearing that threshold)
      -- this is the false-positive rate you pay to keep 95% recall.
    - reject95: threshold = 95th percentile of THIS decoy's scores, so 95% of the decoy is rejected and
      only 5% passes (the "95% precision"/decoy-rejection point). Report the recall of ID-matched
      targets that clears it -- the real-ID recall you keep while holding decoys to a 5% pass-rate.

    Both also carry `precision` = TP/(TP+FP) at that threshold for reference, but the primary numbers are
    recall and decoy_pass, which is what makes the two experiments comparable."""
    if len(pos) == 0 or len(neg) == 0:
        return None, None
    npos, nneg = len(pos), len(neg)

    def stats(t):
        tp = int(np.sum(pos >= t)); fp = int(np.sum(neg >= t))
        return dict(thr=float(t), recall=tp / npos, decoy_pass=fp / nneg,
                    precision=(tp / (tp + fp) if (tp + fp) else 1.0), tp=tp, fp=fp)

    recall95 = stats(float(np.percentile(pos, 5.0)))    # retain 95% of ID-matched targets
    reject95 = stats(float(np.percentile(neg, 95.0)))   # reject 95% of decoys -> 5% pass
    return recall95, reject95


def main():
    summary = {}
    for ds, gt_path, rt_tol, title in DATASETS:
        tgt = load_scores(ds, "target")
        if tgt is None:
            print(f"[{ds}] target scores missing -- skip")
            continue
        gm, gr = load_gt(gt_path)
        matched = match_mask(tgt, gm, gr, rt_tol)
        n_m = int(matched.sum()); n = len(tgt)
        print(f"\n[{ds}] {title}: target n={n}  ID-matched={n_m} ({100*n_m/n:.1f}%)  unmatched={n-n_m}")
        decoy = {d: load_scores(ds, d) for d in DECOYS}
        for d in DECOYS:
            print(f"    {d:9s}: {0 if decoy[d] is None else len(decoy[d])} features")
        summary[ds] = {"title": title, "n_target": n, "n_matched": n_m, "decoys": {}}

        # ---- distributions ----
        fig, axes = plt.subplots(1, 2, figsize=(13, 6.6))
        for ax, (score, xlabel) in zip(axes, SCORES):
            bins = np.linspace(0, 1, 60)
            ax.hist(tgt[score][~matched], bins=bins, density=True, histtype="step", lw=1.8,
                    color=COLORS["unmatched"], label=f"target, unmatched (n={n-n_m})")
            ax.hist(tgt[score][matched], bins=bins, density=True, histtype="stepfilled", lw=2.0, alpha=0.35,
                    edgecolor=COLORS["matched"], facecolor=COLORS["matched"], label=f"target, ID-matched (n={n_m})")
            for d in DECOYS:
                if decoy[d] is None:
                    continue
                _, _, a = roc(tgt[score][matched], decoy[d][score])
                ax.hist(decoy[d][score], bins=bins, density=True, histtype="step", lw=1.8,
                        color=COLORS[d], label=f"{LABELS[d]} (n={len(decoy[d])}, AUC={a:.3f})")
            ax.set_xlabel(xlabel); ax.set_ylabel("density")
            ax.set_title(score.replace("score", "").upper() + " score")
            # legend below the axis (outside the plot) so it never covers the curves
            ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.13),
                      ncol=2, frameon=False)
            ax.spines[["top", "right"]].set_visible(False)
        fig.suptitle(f"{title} - ID-matched target vs decoys", fontsize=13, weight="bold")
        fig.tight_layout(rect=[0, 0.0, 1, 0.95])
        fig.savefig(os.path.join(OUT, f"dist_{ds}.png"), dpi=140, bbox_inches="tight"); plt.close(fig)

        # ---- ROC + AUC + thresholds ----
        figr, axesr = plt.subplots(1, 2, figsize=(11, 5.2))
        for ax, (score, xlabel) in zip(axesr, SCORES):
            pos = tgt[score][matched]
            for d in DECOYS:
                if decoy[d] is None:
                    continue
                neg = decoy[d][score]
                fpr, tpr, a = roc(pos, neg)
                sel = np.linspace(0, len(fpr) - 1, min(len(fpr), 2000)).astype(int)
                ax.plot(fpr[sel], tpr[sel], color=COLORS[d], lw=2.0, label=f"{LABELS[d]}  AUC={a:.3f}")
                rec95, prec95 = op_points(pos, neg)
                node = summary[ds]["decoys"].setdefault(d, {})
                node.setdefault(score, {})["auc"] = a
                node[score]["rec95"] = rec95
                node[score]["prec95"] = prec95
                node["n"] = len(neg)
            ax.plot([0, 1], [0, 1], color="#999999", lw=1.0, ls=":")
            ax.set_xlabel("false-positive rate (decoys accepted)")
            ax.set_ylabel("true-positive rate (matched targets accepted)")
            ax.set_title(xlabel); ax.set_xlim(0, 1); ax.set_ylim(0, 1.001)
            ax.legend(fontsize=8, loc="lower right")
            ax.spines[["top", "right"]].set_visible(False)
        figr.suptitle(f"{title} - ROC: ID-matched target vs decoy", fontsize=13, weight="bold")
        figr.tight_layout(rect=[0, 0, 1, 0.96])
        figr.savefig(os.path.join(OUT, f"roc_{ds}.png"), dpi=140); plt.close(figr)
        print(f"[{ds}] wrote plots/dist_{ds}.png, plots/roc_{ds}.png")

    with open(os.path.join(HERE, "topdown_decoy_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    # ---- text tables ----
    print("\n=== AUC (ID-matched target vs decoy; higher = better separation) ===")
    print(f"{'dataset':<9}{'decoy':<26}{'AUC 2-D':>9}{'AUC 3-D':>9}")
    for ds, s in summary.items():
        for d in DECOYS:
            nd = s["decoys"].get(d)
            if not nd:
                continue
            a2 = nd.get("score2d", {}).get("auc", float("nan"))
            a3 = nd.get("score3d", {}).get("auc", float("nan"))
            print(f"{ds:<9}{LABELS[d]:<26}{a2:>9.3f}{a3:>9.3f}")

    for label, key in [("95% RECALL threshold (5th pct of ID-matched targets: retain 95% of real IDs)", "rec95"),
                       ("95% DECOY-REJECTION threshold (95th pct of decoy: only 5% of decoys pass)", "prec95")]:
        for score, sl in SCORES:
            print(f"\n=== {label} -- {sl.split(' env')[0]} ===")
            print(f"{'dataset':<9}{'decoy':<26}{'thr':>8}{'real recall':>12}{'decoy pass':>12}")
            for ds, s in summary.items():
                for d in DECOYS:
                    nd = s["decoys"].get(d, {}).get(score, {})
                    o = nd.get(key)
                    if not o:
                        continue
                    print(f"{ds:<9}{LABELS[d]:<26}{o['thr']:>8.4f}{100*o['recall']:>11.1f}%{100*o['decoy_pass']:>11.2f}%")
    print("\nwrote topdown_decoy_summary.json")


if __name__ == "__main__":
    main()
