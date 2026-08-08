import json, csv, glob, os
import numpy as np

ci = 3
pi = 0

boltz_predictions = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_Q8/cycle_1/boltz/outputs/boltz_results_input/predictions"
chai_output_path = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_Q8/cycle_1/chai/outputs"
esm3_output_path = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_Q8/cycle_1/esm3/outputs"
af3_output_path = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_Q8/cycle_1/af3/outputs"
trajectory_cycle_dir = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_Q8/design_trajectory/cycle_1"
os.makedirs(trajectory_cycle_dir, exist_ok=True)

def extract_num(s):
    import re
    m = re.search(r'(?:^|_)id(\d+)(?:$|_)', s)
    if m:
        return m.group(1)
    return s

boltz_pl = {}
for seq_dir in glob.glob(f"{boltz_predictions}/*/"):
    sid_raw = os.path.basename(seq_dir.rstrip('/'))
    sid = extract_num(sid_raw)
    vals = []
    for jf in glob.glob(f"{seq_dir}confidence_*.json"):
        with open(jf) as f: d = json.load(f)
        pc = d.get("pair_chains_iptm")
        if pc:
            try: vals.append((pc[str(pi)][str(ci)] + pc[str(ci)][str(pi)])/2)
            except (KeyError, TypeError): pass
    if vals: boltz_pl[sid] = max(vals)

chai_pl = {}
for npz in glob.glob(f"{chai_output_path}/*/scores.model_idx_*.npz"):
    sid_raw = os.path.basename(os.path.dirname(npz))
    sid = extract_num(sid_raw)
    d = np.load(npz)
    if "per_chain_pair_iptm" not in d: continue
    m = d["per_chain_pair_iptm"]
    if m.ndim == 3: m = m[0]
    if ci < m.shape[0] and pi < m.shape[1]:
        s = (m[pi, ci] + m[ci, pi]) / 2
        if sid not in chai_pl or s > chai_pl[sid]: chai_pl[sid] = s

esm3 = {}
ec = f"{esm3_output_path}/esm3_scores.csv"
if os.path.exists(ec):
    with open(ec) as f:
        for row in csv.DictReader(f):
            sid = extract_num(row["id"])
            try: esm3[sid] = float(row["ptm"])
            except (KeyError, ValueError): pass

af3 = {}
files = glob.glob(f"{af3_output_path}/batch_*/*_summary_confidences.json") + \
        glob.glob(f"{af3_output_path}/batch_*/*/*_summary_confidences.json")
for f in files:
    try:
        with open(f) as fh: d = json.load(fh)
        sid = extract_num(os.path.basename(f).split("_summary")[0])
        af3[sid] = d["chain_pair_iptm"][pi][ci]
    except Exception: pass

all_ids = set(chai_pl) | set(af3)
print(f"Boltz ids: {len(boltz_pl)}, ESM3 ids: {len(esm3)}, Chai ids: {len(chai_pl)}, AF3 ids: {len(af3)}, union(chai,af3): {len(all_ids)}")

results = []
for sid in all_ids:
    c = chai_pl.get(sid, 0.0); a = af3.get(sid, 0.0)
    b = boltz_pl.get(sid, 0.0); e = esm3.get(sid, 0.0)
    results.append({"id": sid, "boltz_pl": b, "esm3_ptm": e, "chai_pl": c, "af3_pl": a,
                    "track1_af3only": a, "track2_af3chai": (a+c)/2,
                    "track3_all4": (a+c+b+e)/4})
ar = sorted(results, key=lambda r: -r["af3_pl"]); cr = sorted(results, key=lambda r: -r["chai_pl"])
arank = {r["id"]: i+1 for i,r in enumerate(ar)}; crank = {r["id"]: i+1 for i,r in enumerate(cr)}
for r in results:
    r["track4_rank_sum"] = arank[r["id"]] + crank[r["id"]]
    r["track5_chai_only"] = r["chai_pl"]

out_csv = f"{trajectory_cycle_dir}/all_scores_Q8.csv"
fn = ["id","boltz_pl","esm3_ptm","chai_pl","af3_pl","track1_af3only","track2_af3chai",
      "track3_all4","track4_rank_sum","track5_chai_only"]
with open(out_csv,"w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=fn); w.writeheader(); w.writerows(results)
print(f"Scores for {len(results)} sequences written to {out_csv}")

specs = [("track1_af3only","Track1_AF3",True),("track2_af3chai","Track2_AF3Chai_mean",True),
         ("track3_all4","Track3_All4_mean",True),("track4_rank_sum","Track4_RankSumConsensus",False),
         ("track5_chai_only","Track5_Chai_only",True)]
seeds_csv = f"{trajectory_cycle_dir}/next_cycle_seeds.csv"
seeds_txt = f"{trajectory_cycle_dir}/next_cycle_seeds.txt"
seed_rows = []
lines = [f"RQ2 Q8 cycle 1 -- per-track seed selection", f"From {len(results)} scored sequences", ""]
for k, nm, hib in specs:
    if not results:
        lines.append(f"{nm}: NO RESULTS")
        continue
    best = (max if hib else min)(results, key=lambda r: r[k])
    seed_rows.append({"track": nm, "winning_id": best["id"], "track_score": best[k],
                      "boltz_pl": best["boltz_pl"], "esm3_ptm": best["esm3_ptm"],
                      "chai_pl": best["chai_pl"], "af3_pl": best["af3_pl"]})
    lines.append(f"{nm}: id{best['id']} (score={best[k]:.4f}, boltz={best['boltz_pl']:.3f}, "
                 f"esm3={best['esm3_ptm']:.3f}, chai={best['chai_pl']:.3f}, af3={best['af3_pl']:.3f})")

with open(seeds_csv,"w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=["track","winning_id","track_score","boltz_pl","esm3_ptm","chai_pl","af3_pl"])
    w.writeheader(); w.writerows(seed_rows)

uniq = sorted(set(r["winning_id"] for r in seed_rows), key=lambda x: int(x)) if seed_rows else []
lines.append("")
lines.append(f"Unique seed ids ({len(uniq)}): " + ", ".join(f"id{x}" for x in uniq))
with open(seeds_txt,"w") as f: f.write("\n".join(lines)+"\n")
for l in lines: print(l)
