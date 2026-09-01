# RQ1 C3 pooled compile script - CORRECTED
# dual seed lineage - id40/id141 - Figure 2C

import json, csv, glob, os, re, sys
import numpy as np

WORK = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline"
CYCLE_DIR = os.path.join(WORK, sys.argv[1])
seed_dirs = [os.path.join(WORK, sys.argv[2]), os.path.join(WORK, sys.argv[3])]
TRAJ = os.path.join(WORK, sys.argv[4])
os.makedirs(TRAJ, exist_ok=True)

protein_chain_count = 3
pairs = [(i, j) for i in range(protein_chain_count) for j in range(i + 1, protein_chain_count)]

def extract_orig_id(s):
    m = re.search(r'_id(\d+)$', s)
    return m.group(1) if m else None

all_results = {}

for SEED_DIR in seed_dirs:
    seed_name = os.path.basename(SEED_DIR)
    boltz = {}
    for seq_dir in glob.glob(f"{SEED_DIR}/boltz/outputs/boltz_results_input/predictions/*/"):
        sid = extract_orig_id(os.path.basename(seq_dir.rstrip('/')))
        if not sid: continue
        vals = []
        for jf in glob.glob(f"{seq_dir}confidence_*.json"):
            with open(jf) as f: d = json.load(f)
            pc = d.get("pair_chains_iptm")
            if pc:
                try:
                    v = np.mean([(pc[str(i)][str(j)]+pc[str(j)][str(i)])/2 for i,j in pairs])
                    vals.append(v)
                except: pass
        if vals: boltz[sid] = max(vals)

    esm3 = {}
    esm3_csv = f"{SEED_DIR}/esm3/outputs/esm3_scores.csv"
    if os.path.exists(esm3_csv):
        with open(esm3_csv) as f:
            for row in csv.DictReader(f):
                sid = extract_orig_id(row["id"])
                if sid:
                    try: esm3[sid] = float(row["ptm"])
                    except: pass

    print(f"{seed_name}: Boltz={len(boltz)} ESM3={len(esm3)}")
    for sid in set(boltz)|set(esm3):
        all_results[(seed_name, sid)] = {
            "seed": seed_name, "orig_id": sid,
            "boltz_pp": round(boltz.get(sid,0.0),4),
            "esm3_ptm": round(esm3.get(sid,0.0),4),
            "chai_pp": 0.0, "af3_pp": 0.0,
        }

chai_pattern = re.compile(r'^seq_c5_(?P<seed>id\d+)_id(?P<origid>\d+)_id(?P<pooledidx>\d+)$')
pooled_idx_map = {}
chai_count = 0
for npz in glob.glob(f"{CYCLE_DIR}/chai/outputs/*/scores.model_idx_*.npz"):
    did = os.path.basename(os.path.dirname(npz))
    m = chai_pattern.match(did)
    if not m: continue
    seed, origid, pooledidx = m.group("seed"), m.group("origid"), m.group("pooledidx")
    pooled_idx_map[pooledidx] = (seed, origid)
    d = np.load(npz, allow_pickle=True)
    if "per_chain_pair_iptm" not in d: continue
    mat = d["per_chain_pair_iptm"]
    if mat.ndim == 3: mat = mat[0]
    vals = [(mat[i][j]+mat[j][i])/2 for i,j in pairs if i < mat.shape[0] and j < mat.shape[1]]
    if vals:
        s = float(np.mean(vals))
        key = (seed, origid)
        if key in all_results:
            all_results[key]["chai_pp"] = round(s,4)
        else:
            all_results[key] = {"seed":seed,"orig_id":origid,"boltz_pp":0.0,"esm3_ptm":0.0,"chai_pp":round(s,4),"af3_pp":0.0}
        chai_count += 1
print(f"Chai matched: {chai_count}")

af3_count = 0
for jf in glob.glob(f"{CYCLE_DIR}/af3/outputs/batch_*/id*/*summary_confidences.json"):
    parts = jf.split('/')
    if 'seed-' in parts[-2]: continue
    m = re.search(r'id(\d+)$', parts[-2])
    if not m: continue
    pooledidx = m.group(1)
    if pooledidx not in pooled_idx_map: continue
    seed, origid = pooled_idx_map[pooledidx]
    with open(jf) as f: d = json.load(f)
    pc = d.get("chain_pair_iptm", [])
    if len(pc) >= 3:
        vals = [(pc[i][j]+pc[j][i])/2 for i,j in pairs if i < len(pc) and j < len(pc[i])]
        if vals:
            key = (seed, origid)
            if key in all_results:
                all_results[key]["af3_pp"] = round(float(np.mean(vals)),4)
                af3_count += 1
print(f"AF3 matched: {af3_count}")

results_list = []
for (seed, sid), r in all_results.items():
    b,e,c,a = r["boltz_pp"], r["esm3_ptm"], r["chai_pp"], r["af3_pp"]
    results_list.append({"id":f"{seed}_id{sid}","seed":seed,"orig_id":sid,"boltz_pp":b,"esm3_ptm":e,"chai_pp":c,"af3_pp":a,"combined_all4":round((b+e+c+a)/4,4)})
results_list.sort(key=lambda r: -r["combined_all4"])

print("\nTop 10:")
for r in results_list[:10]:
    print(f"  {r['id']}: combined={r['combined_all4']:.4f} chai={r['chai_pp']:.3f} af3={r['af3_pp']:.3f}")

fn = ["id","seed","orig_id","boltz_pp","esm3_ptm","chai_pp","af3_pp","combined_all4"]
with open(f"{TRAJ}/all_scores_pooled.csv","w",newline="") as f:
    w = csv.DictWriter(f,fieldnames=fn); w.writeheader(); w.writerows(results_list)
print(f"Saved to {TRAJ}/all_scores_pooled.csv")
