#!/usr/bin/env python

# Figure 2D - C3 Trimer: Iterative Design Trajectory

import matplotlib.pyplot as plt
import matplotlib as mpl

cycles = [1, 2, 3, 4, 5]
cycle_labels = ["Cycle 1\n(id14)", "Cycle 2\n(id35)", "Cycle 3\n(id40)", "Cycle 4\n(id121)", "Cycle 5\n(id4)"]

boltz = [0.282, 0.639, 0.356, 0.267, 0.673]
chai = [0.338, 0.452, 0.580, 0.569, 0.298]
af3 = [0.113, 0.247, 0.273, 0.287, 0.123]
esm3 = [0.909, 0.903, 0.927, 0.921, 0.931]

# per-point label offset direction: True = above, False = below
# manually changed based on space, and if say two values were overlapping 
boltz_above = [False, True, True, False, True]
chai_above = [True, False, True, True, False]
af3_above = [False, False, False, True, False]
esm3_above = [True, True, True, True, True]

mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["font.size"] = 13
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.facecolor"] = "white"
mpl.rcParams["figure.facecolor"] = "white"

fig, ax = plt.subplots(figsize=(12, 7.5), dpi=300)
ax.set_facecolor("white")

# colours matching report set ones - i.e. green for af3, orange for chai, etc... 
series = [
    ("Boltz-2 (protein-protein)", boltz, boltz_above, "#4C72B0"),
    ("Chai-1 (protein-protein)", chai, chai_above, "#DD8452"),
    ("AlphaFold 3 (protein-protein)", af3, af3_above, "#55A868"),
    ("ESM3 (pTM)", esm3, esm3_above, "#8172B2"),
]

for label, vals, above_flags, color in series:
    ax.plot(cycles, vals, marker="o", markersize=12, linewidth=2.2,
             color=color, label=label)
    for x, y, above in zip(cycles, vals, above_flags):
        offset = 0.045 if above else -0.045
        va = "bottom" if above else "top"
        ax.text(x, y + offset, f"{y:.3f}", ha="center", va=va,
                 fontsize=11, color=color, fontweight="normal")

ax.set_xticks(cycles)
ax.set_xticklabels(cycle_labels, fontsize=12)
ax.set_xlabel("Design Cycle", fontsize=14)
ax.set_ylabel("ipTM Score", fontsize=14)
ax.set_title("C3 Trimer: Iterative Design Trajectory", fontsize=18, fontweight="bold", pad=20)
ax.set_ylim(0, 1.05)
ax.grid(axis="y", linestyle="--", alpha=0.4, color="#cccccc")
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=True, fontsize=11)

fig.tight_layout()
fig.savefig("rq1_c3_trajectory_fig2d.png", dpi=300, bbox_inches="tight", facecolor="white")
print("Saved: rq1_c3_trajectory_fig2d.png")
