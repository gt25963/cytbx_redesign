#!/usr/bin/env python

# C3 full candidate score distribution per cycle, all 4 tools, single lineage per cycle.
# Cycle 4 = lineage id40. Cycle 5 = lineage id141.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

BASE = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/main_pipeline"
TRUE_IPTM_DIR = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/true_boltz_iptm"

ESM_SOURCES = {
    1: (f"{BASE}/CytbX_4tool_C3/cycle_1/boltz/outputs/combined_scores.csv", "esm3_ptm", None),
    2: (f"{BASE}/CytbX_4tool_C3/cycle_2/boltz/outputs/combined_scores.csv", "esm3_ptm", None),
    3: (f"{BASE}/CytbX_4tool_C3/cycle_3/boltz/outputs/combined_scores.csv", "esm3_ptm", None),
    4: (f"{BASE}/CytbX_4tool_C3_id40/design_trajectory/cycle_4/combined_scores.csv", "esm3_ptm", None),
    5: (f"{BASE}/CytbX_4tool_C3_id40_cycle5/pooled_traj/all_scores_pooled.csv", "esm3_ptm", "cycle_5_seed_id141"),
}

CHAI_AF3_SOURCES = {
    1: (f"{BASE}/CytbX_4tool_C3/cycle_1/c3_cycle1_full_scores.csv", "chai_corrected", "af3"),
    2: (f"{BASE}/CytbX_4tool_C3/design_trajectory/cycle_2/all_scores_c2_compiled.csv", "chai_pp", "af3_pp"),
    3: (f"{BASE}/CytbX_4tool_C3/design_trajectory/cycle_3/all_scores_CORRECTED_full50.csv", "chai_pp", "af3_pp"),
    4: (f"{BASE}/CytbX_4tool_C3_pooled/design_trajectory/cycle_4/all_scores_C3_cycle4_pooled.csv", "chai_pp", "af3_pp"),
}
CHAI_CYCLE5 = f"{BASE}/CytbX_4tool_C3_id40_cycle5/pooled_traj/csvs/chai.csv"
AF3_CYCLE5 = f"{BASE}/CytbX_4tool_C3_id40_cycle5/pooled_traj/csvs/af3.csv"
ID141_TOP50 = f"{BASE}/CytbX_4tool_C3_id40_cycle5/cycle_5_seed_id141/traj/top50.csv"


def load_true_boltz(cycle):
    path = f"{TRUE_IPTM_DIR}/cycle_{cycle}_true_protein_iptm.csv"
    df = pd.read_csv(path)
    return df["true_protein_iptm"].dropna().tolist()


def load_esm(cycle):
    path, esm_col, prefix = ESM_SOURCES[cycle]
    df = pd.read_csv(path)
    if prefix:
        df = df[df["id"].astype(str).str.startswith(prefix)]
    return df[esm_col].dropna().tolist()


def load_chai_af3(cycle):
    if cycle == 5:
        id141_ids = set(pd.read_csv(ID141_TOP50)["id"].astype(str))
        chai_df = pd.read_csv(CHAI_CYCLE5)
        chai_df = chai_df[chai_df["id"].astype(str).isin(id141_ids)]
        af3_df = pd.read_csv(AF3_CYCLE5)
        af3_df = af3_df[af3_df["id"].astype(str).isin(id141_ids)]
        return chai_df["chai_protein_pair_iptm"].dropna().tolist(), af3_df["chain_pair_iptm"].dropna().tolist()
    path, chai_col, af3_col = CHAI_AF3_SOURCES[cycle]
    df = pd.read_csv(path)
    if cycle == 4:
        df = df[(df[chai_col] > 0) & (df[af3_col] > 0)]
    chai_vals = df[chai_col].dropna()
    chai_vals = chai_vals[chai_vals > 0].tolist()
    af3_vals = df[af3_col].dropna()
    af3_vals = af3_vals[af3_vals > 0].tolist()
    return chai_vals, af3_vals


records = []
print("=== Sanity check ===")
for cycle in [1, 2, 3, 4, 5]:
    boltz_vals = load_true_boltz(cycle)
    esm_vals = load_esm(cycle)
    chai_vals, af3_vals = load_chai_af3(cycle)
    print(f"Cycle {cycle}: Boltz(true) n={len(boltz_vals)} mean={sum(boltz_vals)/len(boltz_vals):.3f}, "
          f"ESM3 n={len(esm_vals)}, Chai n={len(chai_vals)}, AF3 n={len(af3_vals)}")
    for v in boltz_vals:
        records.append({"cycle": cycle, "tool": "Boltz-2", "score": v})
    for v in esm_vals:
        records.append({"cycle": cycle, "tool": "ESM3", "score": v})
    for v in chai_vals:
        records.append({"cycle": cycle, "tool": "Chai-1", "score": v})
    for v in af3_vals:
        records.append({"cycle": cycle, "tool": "AlphaFold 3", "score": v})

df_all = pd.DataFrame(records)

mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["font.size"] = 12
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.facecolor"] = "white"
mpl.rcParams["figure.facecolor"] = "white"

TOOL_COLORS = {"Boltz-2": "#4C72B0", "ESM3": "#8172B2", "Chai-1": "#DD8452", "AlphaFold 3": "#55A868"}
TOOLS = ["Boltz-2", "ESM3", "Chai-1", "AlphaFold 3"]
CYCLE_LABELS = {1: "Cycle 1", 2: "Cycle 2", 3: "Cycle 3", 4: "Cycle 4", 5: "Cycle 5"}


def make_plot(kind, outfile):
    # 'box' -> fig 2C
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=300)
    axes = axes.flatten()
    rng = np.random.default_rng(42)
    for ax, tool in zip(axes, TOOLS):
        sub = df_all[df_all.tool == tool]
        cycles = sorted(sub.cycle.unique())
        ymin, ymax = sub["score"].min(), sub["score"].max()
        yrange = ymax - ymin
        pad = max(yrange * 0.35, 0.05)
        ax.set_ylim(max(0, ymin - yrange * 0.1), min(1.05, ymax + pad))

        if kind == "box":
            data = [sub[sub.cycle == c]["score"].values for c in cycles]
            ax.boxplot(data, positions=cycles, widths=0.6, patch_artist=True, showfliers=False,
                       medianprops=dict(color="black", linewidth=1.5),
                       boxprops=dict(facecolor=TOOL_COLORS[tool], alpha=0.6, edgecolor="black", linewidth=0.9),
                       whiskerprops=dict(color="black", linewidth=1.0),
                       capprops=dict(color="black", linewidth=1.0))
            for c in cycles:
                y = sub[sub.cycle == c]["score"].values
                x = rng.normal(c, 0.05, size=len(y))
                ax.scatter(x, y, s=8, color="#333333", alpha=0.25, zorder=3)
        else:
            for c in cycles:
                y = sub[sub.cycle == c]["score"].values
                x = rng.normal(c, 0.06, size=len(y))
                ax.scatter(x, y, s=14, color=TOOL_COLORS[tool], alpha=0.55,
                           edgecolor="black", linewidth=0.3, zorder=3)
                med = np.median(y)
                ax.hlines(med, c - 0.22, c + 0.22, color="black", linewidth=1.5, zorder=4)

        for c in cycles:
            y = sub[sub.cycle == c]["score"].values
            med = np.median(y)
            top = y.max()
            ax.text(c, top + pad * 0.3, f"{med:.2f}", ha="center", va="bottom", fontsize=11,
                     color="#111111", fontweight="bold", zorder=5,
                     bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85))

        ax.set_xticks(cycles)
        ax.set_xticklabels([CYCLE_LABELS[c] for c in cycles], fontsize=10)
        ylabel = "Boltz-2 protein iptm" if tool == "Boltz-2" else f"{tool} score"
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(tool, fontsize=13, fontweight="bold", pad=18)
        ax.grid(axis="y", linestyle="--", alpha=0.4, color="#cccccc")

    fig.suptitle("C3 Trimer: Full Candidate Score Distribution per Cycle, All Tools",
                  fontsize=16, fontweight="bold", y=1.01)
    fig.tight_layout()
    fig.savefig(outfile, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved: {outfile}")


make_plot("strip", "c3_strip_corrected.png") ## not used in report - Paul prefered boxplot
make_plot("box", "c3_box_corrected.png") ## final fig 2C
