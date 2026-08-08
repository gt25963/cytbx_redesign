"""Extract flat id,score CSVs for compile_tracks_v2.py, keyed by pooled index."""
import json, csv, glob, os, re, sys
import numpy as np

WORK = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline"
CYCLE_DIR = os.path.join(WORK, sys.argv[1])
seed_dirs = [os.path.join(WORK, sys.argv[2]), os.path.join(WORK, sys.argv[3])]
OUT = os.path.join(WORK, sys.argv[4])
os.makedirs(OUT, exist_ok=True)

pairs = [(0,1),(0,2),(1,2)]

def orig_id(s):
    m = re.search(r'_id(\d+)$', s)
    return m.group(1) if m else None

chai_pattern = re.compile(r'^seq_c5_(?P<seed>id\d+)_id(?P<origid>\d+)_id(?P<pooledidx>\d+)$')
pooled_idx_map = {}
chai_best = {}
for npz in glob.glob(f"{CYCLE_DIR}/chai/outputs/*/scores.model_idx_*.npz"):
    did = os.path.basename(os.path.dirname(npz))
    m = chai_pattern.match(did)
    if not m: continue
    seed, oid, pidx = m.group("seed"), m.group("origid"), m.group("pooledidx")
    pooled_idx_map[(seed, oid)] = pidx
    d = np.load(npz, allow_pickle=True)
    if "per_chain_pair_iptm" not in d: continue
    mat = d["per_chain_pair_iptm"]
    if mat.ndim == 3: mat = mat[0]
    vals = [(mat[i][j]+mat[j][i])/2 for i,j in pairs if i < mat.shape[0] and j < mat.shape[1]]
    if vals:
        s = round(float(np.mean(vals)),4)
        if pidx not in chai_best or s > chai_best[pidx]:
            chai_best[pidx] = s
chai_rows = [{"id": pidx, "chai_protein_pair_iptm": s} for pidx, s in chai_best.items()]
print(f"Chai: {len(chai_rows)} rows, {len(pooled_idx_map)} pooled-idx mappings")

with open(f"{OUT}/chai.csv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=["id","chai_protein_pair_iptm"]); w.writeheader(); w.writerows(chai_rows)

af3_rows = []
for jf in glob.glob(f"{CYCLE_DIR}/af3/outputs/batch_*/id*/*summary_confidences.json"):
    parts = jf.split('/')
    if 'seed-' in parts[-2]: continue
    m = re.search(r'id(\d+)$', parts[-2])
    if not m: continue
    pidx = m.group(1)
    with open(jf) as f: d = json.load(f)
    pc = d.get("chain_pair_iptm", [])
    if len(pc) >= 3:
        vals = [(pc[i][j]+pc[j][i])/2 for i,j in pairs if i < len(pc) and j < len(pc[i])]
        if vals:
            af3_rows.append({"id": pidx, "chain_pair_iptm": round(float(np.mean(vals)),4)})
print(f"AF3: {len(af3_rows)} rows")

with open(f"{OUT}/af3.csv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=["id","chain_pair_iptm"]); w.writeheader(); w.writerows(af3_rows)

boltz_rows, esm_rows = [], []
for SEED_DIR in seed_dirs:
    full_name = os.path.basename(SEED_DIR)
    m_seed = re.search(r'(id\d+)$', full_name)
    seed_name = m_seed.group(1) if m_seed else full_name
    for seq_dir in glob.glob(f"{SEED_DIR}/boltz/outputs/boltz_results_input/predictions/*/"):
        oid = orig_id(os.path.basename(seq_dir.rstrip('/')))
        if not oid: continue
        pidx = pooled_idx_map.get((seed_name, oid))
        if not pidx: continue
        vals = []
        for jf in glob.glob(f"{seq_dir}confidence_*.json"):
            with open(jf) as f: d = json.load(f)
            pc = d.get("pair_chains_iptm")
            if pc:
                try:
                    vals.append(np.mean([(pc[str(i)][str(j)]+pc[str(j)][str(i)])/2 for i,j in pairs]))
                except: pass
        if vals:
            boltz_rows.append({"id": pidx, "iptm": round(float(max(vals)),4)})

    esm3_csv = f"{SEED_DIR}/esm3/outputs/esm3_scores.csv"
    if os.path.exists(esm3_csv):
        with open(esm3_csv) as f:
            for row in csv.DictReader(f):
                oid = orig_id(row["id"])
                if not oid: continue
                pidx = pooled_idx_map.get((seed_name, oid))
                if not pidx: continue
                try:
                    esm_rows.append({"id": pidx, "ptm": round(float(row["ptm"]),4)})
                except: pass
print(f"Boltz: {len(boltz_rows)} rows, ESM3: {len(esm_rows)} rows")

with open(f"{OUT}/boltz.csv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=["id","iptm"]); w.writeheader(); w.writerows(boltz_rows)
with open(f"{OUT}/esm.csv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=["id","ptm"]); w.writeheader(); w.writerows(esm_rows)

print(f"All four CSVs written to {OUT}")
