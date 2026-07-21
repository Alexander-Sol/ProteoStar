#!/usr/bin/env python3
"""Visual example of the target averagine and the three decoy envelopes at a top-down-relevant mass.

Reads envelopes_12kDa.tsv (model, k, weight) dumped by the `dump_envelopes` Rust example, which uses
EnvelopeModel::comb_weights keyed by the most-intense (mode) mass and normalized to max 1.0. Every model
is therefore drawn with its **most-abundant peak at the same mass (12,000 Da) and height 1** -- only the
envelope shape differs. The `shifted` decoy is not a distinct envelope model (it is the real averagine
weights on a 0.94-Da lattice), so it is derived here from the averagine row with compressed tooth spacing.

Usage:  python plot_averagine_models.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
MODE_MASS = 12000.0
C13 = 1.0033548
SPACING_SCALE = 0.9368  # shifted-spacing decoy: 0.94 Da teeth (= C13 * 0.9368)

COLORS = {"averagine": "#111111", "shifted": "#0072B2", "weird": "#D55E00", "shuffled": "#009E73"}
TITLES = {
    "averagine": "target averagine (real peptide/proteoform envelope)",
    "shifted":   "shifted spacing (0.94 Da) — averagine weights, off-lattice",
    "weird":     "weird averagine (Fe+Cl) — exotic composition, A+2 ladder",
    "shuffled":  "shuffled envelope — averagine weights randomly permuted",
}
ORDER = ["averagine", "shifted", "weird", "shuffled"]


def load():
    rows = np.genfromtxt(os.path.join(HERE, "envelopes_12kDa.tsv"), delimiter="\t",
                         names=True, dtype=None, encoding="utf-8")
    out = {}
    for name in ("averagine", "weird", "shuffled"):
        sel = rows["model"] == name
        k = rows["k"][sel].astype(int)
        w = rows["weight"][sel].astype(float)
        order = np.argsort(k)
        out[name] = w[order] / w.max()
    return out


def teeth(weights, spacing):
    """Place teeth on a neutral-mass axis with the mode (tallest tooth) at MODE_MASS."""
    mode = int(np.argmax(weights))
    x = MODE_MASS + (np.arange(len(weights)) - mode) * spacing
    return x, weights


def main():
    w = load()
    env = {
        "averagine": teeth(w["averagine"], C13),
        "shifted":   teeth(w["averagine"], C13 * SPACING_SCALE),   # same weights, compressed lattice
        "weird":     teeth(w["weird"], C13),
        "shuffled":  teeth(w["shuffled"], C13),
    }
    # shared x-window covering the widest envelope (weird), padded
    xmin = min(x.min() for x, _ in env.values()) - 1.5
    xmax = max(x.max() for x, _ in env.values()) + 1.5

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 7.2), sharex=True, sharey=True)
    for ax, name in zip(axes.ravel(), ORDER):
        x, y = env[name]
        ax.axvline(MODE_MASS, color="#cccccc", lw=1.0, ls="--", zorder=0)
        ml, sl, bl = ax.stem(x, y, basefmt=" ")
        plt.setp(sl, color=COLORS[name], lw=1.6)
        plt.setp(ml, color=COLORS[name], markersize=3.5)
        ax.set_title(TITLES[name], fontsize=10, color=COLORS[name], weight="bold")
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(0, 1.08)
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(0.015, 0.93, f"{len(y)} teeth", transform=ax.transAxes, fontsize=8, color="#555555")
    for ax in axes[-1]:
        ax.set_xlabel("neutral mass (Da)")
    for ax in axes[:, 0]:
        ax.set_ylabel("relative intensity")
    fig.suptitle("Isotope envelope of each model at 12,000 Da\n(most-abundant peak aligned to 12,000 Da and normalized to 1; dashed line = common mode)",
                 fontsize=12.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(HERE, "plots", "averagine_models.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
