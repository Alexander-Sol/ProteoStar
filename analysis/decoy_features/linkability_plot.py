#!/usr/bin/env python3
"""Grouped-bar figure of MS2-linkability rates from ms2_linkability.json.

Three criteria (x, weak->strong): window overlaps a charge state / isolation centre on any isotope /
isolation centre on the most-abundant isotope. Three series per file: ID-matched targets (real, positive
control), unmatched targets (the question), weird-averagine decoys (coincidence null)."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
rows = json.load(open(os.path.join(HERE, "ms2_linkability.json")))
FILES = [("10", "10-min (CA/Lumos)"), ("65", "65-min (glyco)"), ("2h", "2-hr (IonStar)")]
CRIT = [("linkable_pct", "window overlaps\na charge state"),
        ("targeted_pct", "centre on\nany isotope"),
        ("center_pct", "centre on most-\nabundant isotope")]
SERIES = [("matched_target", "ID-matched target (real)", "#111111"),
          ("unmatched_target", "unmatched target (?)", "#D55E00"),
          ("weird_decoy", "weird decoy (null)", "#56B4E9")]

def get(fk, sname):
    for r in rows:
        if r["file"] == fk and r["set"] == sname:
            return r
    return None

fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6), sharey=True)
x = np.arange(len(CRIT)); w = 0.26
for ax, (fk, title) in zip(axes, FILES):
    for j, (sname, slabel, color) in enumerate(SERIES):
        r = get(fk, sname)
        vals = [r[c] for c, _ in CRIT]
        bars = ax.bar(x + (j - 1) * w, vals, w, color=color, label=slabel if fk == "10" else None)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 1.0, f"{v:.0f}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels([lab for _, lab in CRIT], fontsize=8)
    ax.set_title(title, fontsize=11, weight="bold")
    ax.set_ylim(0, 100)
    ax.spines[["top", "right"]].set_visible(False)
    if fk == "10":
        ax.set_ylabel("% of features above threshold")
axes[0].legend(fontsize=8, loc="upper right", frameon=False)
fig.suptitle("MS2-linkability of above-threshold features (95%-decoy-reject cut), by criterion",
             fontsize=12.5, weight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])
p = os.path.join(HERE, "plots", "ms2_linkability.png")
fig.savefig(p, dpi=140)
print("wrote", p)
