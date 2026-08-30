#!/usr/bin/env python3

# Simple rescoring using mean Boltz confidence and mean ESM3 pTM per structure.
# No ID matching required.
# Dense-pose version of rescore_simple.py.
## Same mean-based method, applied to denser RPXDock sampling control (oligomer_screen_dense).
## Produces the dense comparison values reported in Results 1.1: 0.72 C2, 0.65 C3 vs the default density baseline (0.76/0.67)

import csv
import os
import glob

screen_dir = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/oligomer_screen_dense"
output_file = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/screen_scores_final_dense.csv"

results = []

for state in [2, 3, 4, 5]:
    ## dense variant only has one pose per state (rpxdock_top_n=1) -> j only ranges over [1] here
    for j in [1]:
        boltz_pred = f"{screen_dir}/C{state}/top{j}/boltz/outputs/boltz_results_input/predictions"
        esm3_out = f"{screen_dir}/C{state}/top{j}/esm3/outputs/esm3_scores.csv"

        # mean boltz confidence across every sequence x model in the batch
        boltz_scores = []
        for seq_dir in glob.glob(f"{boltz_pred}/*/"):
            for json_file in glob.glob(f"{seq_dir}confidence_*.json"):
                import json
                with open(json_file) as f:
                    data = json.load(f)
                boltz_scores.append(data.get("confidence_score", 0.0))
        mean_boltz = sum(boltz_scores) / len(boltz_scores) if boltz_scores else 0.0

        # mean esm3 ptm across every sequence in the batch
        esm3_ptms = []
        if os.path.exists(esm3_out):
            with open(esm3_out) as f:
                for row in csv.DictReader(f):
                    try:
                        esm3_ptms.append(float(row["ptm"]))
                    except:
                        pass
        mean_esm3 = sum(esm3_ptms) / len(esm3_ptms) if esm3_ptms else 0.0

        combined = (mean_boltz + mean_esm3) / 2.0
        results.append((f"C{state}", f"top{j}", mean_boltz, mean_esm3, combined))
        print(f"C{state} top{j}: boltz={mean_boltz:.4f} esm3={mean_esm3:.4f} combined={combined:.4f}")

results.sort(key=lambda x: x[4], reverse=True)

with open(output_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["state", "top_n", "mean_boltz", "mean_esm3", "combined"])
    for r in results:
        writer.writerow(r)

print(f"\nBest: {results[0][0]} {results[0][1]} with combined score {results[0][4]:.4f}")
print(f"Scores saved to {output_file}")
