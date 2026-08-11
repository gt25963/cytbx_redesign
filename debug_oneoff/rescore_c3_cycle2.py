import json, csv, glob, os, re
import numpy as np

#compile corrected per-track scores for C3 cycle 2, using true pairwise-averaged boltz confidence extracted directly from raw JSONs, rather than the inflated top-level iptm field (see Methods: "For Boltz-2 the iptm field was used" note and the true-pairwise-Boltz correction described in discussion of report)
base = "main_pipeline/CytbX_4tool_C3/cycle_2"
traj = "main_pipeline/CytbX_4tool_C3/design_trajectory/cycle_2"
pairs = [(0,1),(0,2),(1,2)] # C3 trimer protein pairs

def keyid(s):
    #pull the numeric design id out of any filename/dirname
    m = re.search(r'id(\d+)', s)
    return m.group(1) if m else s

#Boltz: average all three pairwise interface scores per model, keep best model per candidate
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

#Chai: same pairwise averaging, read from per_chain_pair_iptm in the raw npz (not the summary aggregate_score - Methods)
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

#AF3: average the same three pairwise scores from chain_pair_iptm
af3 = {}
for f in glob.glob(f"{base}/af3/outputs/batch_*/*_summary_confidences.json") + glob.glob(f"{base}/af3/outputs/batch_*/*/*_summary_confidences.json"):
    try:
        d = json.load(open(f)); sid = keyid(os.path.basename(f))
        v = [ d["chain_pair_iptm"][i][j] for i,j in pairs ]
        af3[sid] = sum(v)/len(v)
    except Exception: continue

#ESM3: pull pTM directly from the compiled scores CSV for this cycle
esm = {}
ec = f"{base}/esm3/outputs/esm3_scores.csv"
if os.path.exists(ec):
    for row in csv.DictReader(open(ec)):
        sid = keyid(row["id"])
        try: esm[sid] = float(row["ptm"])
        except: pass

#Combine all four tools' scores into one row per candidate, computing the five selection "tracks" used for seed advancement each cycle
all_ids = set(chai) | set(af3)
rows=[]
for sid in all_ids:
    c=chai.get(sid,0.0); a=af3.get(sid,0.0); b=boltz.get(sid,0.0); e=esm.get(sid,0.0)
    rows.append({"id":sid,"boltz_pp":b,"esm3_ptm":e,"chai_pp":c,"af3_pp":a,
                 "track1_af3only":a,"track2_af3chai":(a+c)/2,"track3_all4":(a+c+b+e)/4,
                 "track5_chai_only":c})

#Track4 = combined rank-sum of AF3 and Chai-1 rankings (lower is better)
ar=sorted(rows,key=lambda r:-r["af3_pp"]); cr=sorted(rows,key=lambda r:-r["chai_pp"])
arank={r["id"]:i+1 for i,r in enumerate(ar)}; crank={r["id"]:i+1 for i,r in enumerate(cr)}
for r in rows: r["track4_rank_sum"]=arank[r["id"]]+crank[r["id"]]

print(f"rescored {len(rows)} designs; boltz nonzero: {sum(1 for r in rows if r['boltz_pp']>0)}")

#report the winning candidate under each of the five selection tracks
specs=[("track1_af3only","Track1_AF3",True),("track2_af3chai","Track2_AF3Chai",True),
       ("track3_all4","Track3_All4",True),("track4_rank_sum","Track4_RankSum",False),
       ("track5_chai_only","Track5_Chai",True)]
print("\n=== corrected seeds ===")
for k,name,hi in specs:
    best=(max if hi else min)(rows,key=lambda r:r[k])
    print(f"{name}: id{best['id']}  (track={best[k]:.4f} boltz={best['boltz_pp']:.3f} esm3={best['esm3_ptm']:.3f} chai={best['chai_pp']:.3f} af3={best['af3_pp']:.3f})")

#write the full corrected per-candidate scores table for this cycle
with open(f"{traj}/all_scores_CORRECTED_boltzfix.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["id","boltz_pp","esm3_ptm","chai_pp","af3_pp","track1_af3only","track2_af3chai","track3_all4","track4_rank_sum","track5_chai_only"])
    w.writeheader(); w.writerows(rows)
print(f"\nwrote {traj}/all_scores_CORRECTED_boltzfix.csv")
