#!/usr/bin/env python3
import json, glob, os, csv

base = "rq2/master_pipeline/RQ2_Q8/cycle_1/boltz/outputs/boltz_results_input/predictions"

rows = []
for sd in sorted(glob.glob(f"{base}/*/")):
    sid = os.path.basename(sd.rstrip('/'))
    model_scores = []
    for jf in sorted(glob.glob(f"{sd}confidence_*.json")):
        d = json.load(open(jf))
        pc = d.get("pair_chains_iptm")
        if not pc:
            continue
        v = (pc["0"]["3"] + pc["3"]["0"]) / 2.0
        model_scores.append(v)
    if model_scores:
        rows.append((sid, len(model_scores), round(min(model_scores), 4), round(max(model_scores), 4)))

with open("rq2/master_pipeline/RQ2_Q8/cycle_1/boltz/outputs/q8_cycle1_boltz_scores.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["id", "n_models", "model_min", "boltz_score"])
    for r in rows:
        w.writerow(r)

print(f"wrote {len(rows)} rows")
