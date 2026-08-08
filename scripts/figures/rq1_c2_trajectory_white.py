#!/usr/bin/env python
"""
RQ1 C2 Dimer: Iterative Design Trajectory, white background, 4 cycles.
Values from Martina's corrected scores table.
"""

import matplotlib.pyplot as plt
import matplotlib as mpl

cycles = [1, 2, 3, 4]
cycle_labels = ["Cycle 1\n(id96)", "Cycle 2\n(id66)", "Cycle 3\n(id177)", "Cycle 4\n(id122)"]

boltz = [0.423, 0.531, 0.358, 0.582]
chai = [0.577, 0.275, 0.187, 0.129]
af3 = [0.140, 0.150, 0.200, 0.250]
esm3 = [0.926, 0.938, 0.932, 0.895]

boltz_above = [False, True, False, True]
chai_above = [True, True, False, False]
af3_above = [False, False, True, True]
esm3_above = [True, True, True, True]

mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["font.size"] = 13
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.facecolor"] = "white"
mpl.rcParams["figure.facecolor"] = "white"

fig, ax = plt.subplots(figsize=(11, 7.5), dpi=300)
ax.set_facecolor("white")

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
ax.set_title("RQ1 C2 Dimer: Iterative Design Trajectory", fontsize=18, fontweight="bold", pad=20)
ax.set_ylim(0, 1.05)
ax.grid(axis="y", linestyle="--", alpha=0.4, color="#cccccc")
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=True, fontsize=11)

fig.tight_layout()
fig.savefig("rq1_c2_trajectory_white.png", dpi=300, bbox_inches="tight", facecolor="white")
print("Saved: rq1_c2_trajectory_white.png")
