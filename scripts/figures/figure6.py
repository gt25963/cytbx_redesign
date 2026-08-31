#!/usr/bin/env python

# Figure 6 - panels A - D (C2, C3, FMN, Q8).
# Scatter of every candidate per cycle, Chai-1 vs AF3, with mean values.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv

BASE = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/main_pipeline"
RQ2_BASE = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline"
GAP_DIR = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/gap_cycles"

def load_pairs(path, chai_col, af3_col, id_col="id", prefix_filter=None):
    df = pd.read_csv(path)
    if prefix_filter:
        df = df[df[id_col].astype(str).str.startswith(prefix_filter)]
    df = df[(df[chai_col] > 0) & (df[af3_col] > 0)]
    return df[chai_col].tolist(), df[af3_col].tolist()

def load_gap(name):
    chai_v, af3_v = [], []
    with open(f"{GAP_DIR}/{name}.csv") as f:
        for row in csv.DictReader(f):
            chai_v.append(float(row["chai"]))
            af3_v.append(float(row["af3"]))
    return chai_v, af3_v

c2_cycles = {
    1: load_pairs(f"{BASE}/CytbX_4tool/cycle_1/c2_cycle1_full_scores.csv", "chai", "af3"),
    2: load_pairs(f"{BASE}/CytbX_4tool/design_trajectory/cycle_2/all_scores_c2_compiled.csv", "chai_pp", "af3_pp"),
    3: load_gap("c2_cycle3"),
    4: load_pairs(f"{BASE}/CytbX_4tool_id177/design_trajectory/cycle_4/all_scores_CORRECTED.csv", "chai_pp", "af3_pp"),
}

c3_cycles = {
    1: load_pairs(f"{BASE}/CytbX_4tool_C3/cycle_1/c3_cycle1_full_scores.csv", "chai_corrected", "af3"),
    2: load_pairs(f"{BASE}/CytbX_4tool_C3/design_trajectory/cycle_2/all_scores_c2_compiled.csv", "chai_pp", "af3_pp"),
    3: load_pairs(f"{BASE}/CytbX_4tool_C3/design_trajectory/cycle_3/all_scores_CORRECTED_full50.csv", "chai_pp", "af3_pp"),
    4: load_gap("c3_cycle4"),
    5: load_gap("c3_cycle5"),
}

fmn_cycles = {
    1: load_gap("fmn_cycle1"),
    2: load_pairs(f"{RQ2_BASE}/RQ2_FMN_hemB/design_trajectory/cycle_2/all_scores_FMN_hemB_cycle2.csv", "chai_pl", "af3_pl"),
    3: load_pairs(f"{RQ2_BASE}/RQ2_FMN/design_trajectory/cycle_3/all_scores_FMN_cycle3.csv", "chai_pl", "af3_pl"),
    4: load_pairs(f"{RQ2_BASE}/RQ2_FMN/design_trajectory/cycle_4/all_scores_FMN_cycle4.csv", "chai_pl", "af3_pl"),
    5: load_pairs(f"{RQ2_BASE}/RQ2_FMN/design_trajectory/cycle_5/all_scores_FMN_cycle5_REAL.csv", "chai_pl", "af3_pl"),
}

q8_cycles = {
    1: load_gap("q8_cycle1"),
    2: load_pairs(f"{RQ2_BASE}/RQ2_Q8/design_trajectory/cycle_2/all_scores_Q8_cycle2.csv", "chai_pl", "af3_pl"),
    3: load_pairs(f"{RQ2_BASE}/RQ2_Q8/design_trajectory/cycle_3/all_scores_Q8_cycle3.csv", "chai_pl", "af3_pl"),
    4: load_pairs(f"{RQ2_BASE}/RQ2_Q8/design_trajectory/cycle_4/all_scores_Q8_cycle4_REAL.csv", "chai_pl", "af3_pl"),
    5: load_pairs(f"{RQ2_BASE}/RQ2_Q8/design_trajectory/cycle_5/all_scores_Q8_cycle5_full.csv", "chai_pl", "af3_pl"),
}

mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["font.size"] = 13
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.facecolor"] = "white"
mpl.rcParams["figure.facecolor"] = "white"

def make_track_figure(title, cycles, chai_label, af3_label, outfile):
    rng = np.random.default_rng(42)
    cyc_nums = sorted(cycles.keys())

    fig, ax = plt.subplots(figsize=(9.5, 7), dpi=300)

    all_vals = [v for c in cycles.values() for pair in c for v in pair]
    ymax = max(all_vals)
    top_pad = max(ymax * 0.08, 0.03)
    ax.set_ylim(0, 0.75)

    label_y = ymax * 0.035 + 0.015  # small fixed height just above y=0

    for cyc in cyc_nums:
        chai_v, af3_v = cycles[cyc]
        x_chai = rng.normal(cyc - 0.12, 0.05, size=len(chai_v))
        x_af3 = rng.normal(cyc + 0.12, 0.05, size=len(af3_v))
        ax.scatter(x_chai, chai_v, s=22, color="#DD8452", alpha=0.55,
                   edgecolor="black", linewidth=0.3, zorder=3,
                   label=chai_label if cyc == cyc_nums[0] else None)
        ax.scatter(x_af3, af3_v, s=22, color="#55A868", alpha=0.55,
                   edgecolor="black", linewidth=0.3, zorder=3,
                   label=af3_label if cyc == cyc_nums[0] else None)

        chai_mean = np.mean(chai_v)
        af3_mean = np.mean(af3_v)
        ax.hlines(chai_mean, cyc - 0.22, cyc - 0.02, color="#8B4A2B", linewidth=2, zorder=4)
        ax.hlines(af3_mean, cyc + 0.02, cyc + 0.22, color="#2F5233", linewidth=2, zorder=4)

        ax.text(cyc - 0.12, label_y, f"{chai_mean:.2f}",
                ha="center", va="bottom", fontsize=10, color="#DD8452", fontweight="bold")
        ax.text(cyc + 0.12, label_y, f"{af3_mean:.2f}",
                ha="center", va="bottom", fontsize=10, color="#55A868", fontweight="bold")

    ax.set_xticks(cyc_nums)
    ax.set_xlabel("Design Cycle", fontsize=13)
    ax.set_ylabel("ipTM Score", fontsize=13)
    ax.set_title(title, fontsize=15, fontweight="bold", pad=15)
    ax.grid(axis="y", linestyle="--", alpha=0.4, color="#cccccc")
    ax.legend(loc="upper left", frameon=True, fontsize=11)

    fig.tight_layout()
    fig.savefig(outfile, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved: {outfile}")

make_track_figure("C2 Dimer: Chai-1 vs AlphaFold 3", c2_cycles,
                   "Chai-1 protein-protein ipTM", "AlphaFold 3 protein-protein ipTM",
                   "fig5_c2_scatter_v9.png")
make_track_figure("C3 Trimer: Chai-1 vs AlphaFold 3", c3_cycles,
                   "Chai-1 protein-protein ipTM", "AlphaFold 3 protein-protein ipTM",
                   "fig5_c3_scatter_v9.png")
make_track_figure("FMN: Chai-1 vs AlphaFold 3", fmn_cycles,
                   "Chai-1 protein-cofactor ipTM", "AlphaFold 3 protein-cofactor ipTM",
                   "fig5_fmn_scatter_v9.png")
make_track_figure("Q8: Chai-1 vs AlphaFold 3", q8_cycles,
                   "Chai-1 protein-cofactor ipTM", "AlphaFold 3 protein-cofactor ipTM",
                   "fig5_q8_scatter_v9.png")
