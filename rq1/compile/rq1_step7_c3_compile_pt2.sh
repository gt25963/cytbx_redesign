#!/usr/bin/env python3

# Compiles Boltz-2/ESM3/Chai-1/AF3 scores for every design in a cycle and reports the top design under 5 different track definitions. 
# Called from rq1_step7_c3_compile_pt1.sh with base/traj passed in as arguments.

import json, csv, glob, os, re, sys
import numpy as np

base = sys.argv[1]
traj = sys.argv[2]

# C3 trimer: all three pairwise chain combinations count toward the protein-protein interface score (unlike RQ2, which only has one protein-cofactor pair)
pairs = [(0,1), (0,2), (1,2)]

def keyid(s):
    m = re.search(r'id(\d+)', s)
    return m.group(1) if m else s

# Boltz: average the three pairwise ipTMs per model, then take the best model per design
boltz = {}
for sd in glob.glob(f"{base}/boltz/outputs/boltz_results_input/predictions/*/"):
    sid = keyid(os.path.basename(sd.rstrip('/')))
    vals_models = []
    for jf in glob.glob(f"{sd}confidence_*.json"):
        d = json.load(open(jf)); pc = d.get("pair_chains_iptm")
        if not pc: continue
        v = [ (pc[str(i)][str(j)]+pc[str(j)][str(i)])/2 for i,j in pairs ]
        vals_models.append(sum(v)/len(v))
    if vals_models: boltz[sid] = max(vals_models)

# Chai: same three-pair average, read directly from the per_chain_pair_iptm npz matrix rather than any summary/aggregate score
chai = {}
for npz in glob.glob(f"{base}/chai/outputs/*/scores.model_idx_*.npz"):
    sid = keyid(os.path.basename(os.path.dirname(npz)))
    d = np.load(npz)
    if "per_chain_pair_iptm" not in d: continue
    m = d["per_chain_pair_iptm"]
    if m.ndim == 3: m = m[0]
    v = [ (m[i,j]+m[j,i])/2 for i,j in pairs if i<m.shape[0] and j<m.shape[1] ]
    if not v: continue
    sc = sum(v)/len(v)
    if sid not in chai or sc > chai[sid]: chai[sid] = sc

# AF3: two possible output layouts are globbed (flat batch folder vs per-id subfolder), since AF3's output structure varied across runs
af3 = {}
for f in glob.glob(f"{base}/af3/outputs/batch_*/*_summary_confidences.json") + glob.glob(f"{base}/af3/outputs/batch_*/*/*_summary_confidences.json"):
    try:
        d = json.load(open(f)); sid = keyid(os.path.basename(f))
        v = [ d["chain_pair_iptm"][i][j] for i,j in pairs ]
        af3[sid] = sum(v)/len(v)
    except Exception: continue

# ESM3: single whole-fold pTM per design (no per-chain-pair breakdown available)
esm = {}
ec = f"{base}/esm3/outputs/esm3_scores.csv"
if os.path.exists(ec):
    for row in csv.DictReader(open(ec)):
        sid = keyid(row["id"])
        try: esm[sid] = float(row["ptm"])
        except: pass

# Union of Chai/AF3 ids as the base set - Boltz/ESM3 are fast filters that should already cover every design, but Chai/AF3 define which ones actually got scored
all_ids = set(chai) | set(af3)
rows = []
for sid in all_ids:
    c = chai.get(sid, 0.0); a = af3.get(sid, 0.0); b = boltz.get(sid, 0.0); e = esm.get(sid, 0.0)
    rows.append({"id": sid, "boltz_pp": b, "esm3_ptm": e, "chai_pp": c, "af3_pp": a,
                 "track1_af3only": a, "track2_af3chai": (a+c)/2, "track3_all4": (a+c+b+e)/4,
                 "track5_chai_only": c})

# Track4 (rank-sum) needs both rankings computed before it can be added to each row
ar = sorted(rows, key=lambda r: -r["af3_pp"]); cr = sorted(rows, key=lambda r: -r["chai_pp"])
arank = {r["id"]: i+1 for i,r in enumerate(ar)}; crank = {r["id"]: i+1 for i,r in enumerate(cr)}
for r in rows: r["track4_rank_sum"] = arank[r["id"]] + crank[r["id"]]

print(f"rescored {len(rows)} designs; boltz nonzero: {sum(1 for r in rows if r['boltz_pp']>0)}")

# Five track definitions matching the FMN/Q8 compile scripts' pattern - rank_sum is the only "lower is better" metric, hence the hi=False flag
specs = [("track1_af3only","Track1_AF3",True), ("track2_af3chai","Track2_AF3Chai",True),
         ("track3_all4","Track3_All4",True), ("track4_rank_sum","Track4_RankSum",False),
         ("track5_chai_only","Track5_Chai",True)]

print("\n=== RQ1 cycle_2 seed selection ===")
for k,name,hi in specs:
    best = (max if hi else min)(rows, key=lambda r: r[k])
    print(f"{name}: id{best['id']}  (track={best[k]:.4f} boltz={best['boltz_pp']:.3f} esm3={best['esm3_ptm']:.3f} chai={best['chai_pp']:.3f} af3={best['af3_pp']:.3f})")

with open(f"{traj}/all_scores_c2_compiled.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["id","boltz_pp","esm3_ptm","chai_pp","af3_pp","track1_af3only","track2_af3chai","track3_all4","track4_rank_sum","track5_chai_only"])
    w.writeheader(); w.writerows(rows)
print(f"\nwrote {traj}/all_scores_c2_compiled.csv")
