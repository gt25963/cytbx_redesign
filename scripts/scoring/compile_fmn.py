"""
RQ2 FMN compile script.
Usage: python compile_fmn.py <cycle_number>
Example: python compile_fmn.py 5
"""
import json, csv, glob, os, re, sys
import numpy as np
from Bio.PDB import MMCIFParser

WORK = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline"
RQ2_FMN = f"{WORK}/rq2/master_pipeline/RQ2_FMN"
cycle = sys.argv[1]
CYCLE_DIR = f"{RQ2_FMN}/cycle_{cycle}"
TRAJ = f"{RQ2_FMN}/design_trajectory/cycle_{cycle}"
os.makedirs(TRAJ, exist_ok=True)

pi = 0
ci_fmn = 2

def extract_orig_id(s):
    m = re.search(r'_id(\d+)$', s)
    return m.group(1) if m else None

# Load top50 mapping
top50 = []
with open(f"{TRAJ}/top50_FMN_cycle{cycle}.csv") as f:
    for row in csv.DictReader(f):
        top50.append((row["lineage"], row["original_id"]))
id_mapping = {str(i+1): top50[i] for i in range(len(top50))}

seeds = list(set(s for s,_ in top50))
print(f"Seeds: {seeds}")

# Boltz
print("Extracting Boltz...")
boltz_pl = {}
for seed in seeds:
    pred_dir = f"{CYCLE_DIR}/seed_{seed}/boltz_outputs/boltz_results_input/predictions"
    for seq_dir in glob.glob(f"{pred_dir}/*/"):
        orig_id = extract_orig_id(os.path.basename(seq_dir.rstrip('/')))
        if not orig_id: continue
        vals = []
        for jf in glob.glob(f"{seq_dir}confidence_*.json"):
            with open(jf) as f: d = json.load(f)
            pc = d.get("pair_chains_iptm")
            if pc:
                try: vals.append((pc[str(pi)][str(ci_fmn)]+pc[str(ci_fmn)][str(pi)])/2)
                except: pass
        if vals:
            for nid,(s,oid) in id_mapping.items():
                if s==seed and oid==orig_id: boltz_pl[nid]=max(vals)
print(f"  Boltz: {len(boltz_pl)}")

# ESM3
print("Extracting ESM3...")
esm3 = {}
for seed in seeds:
    esm3_csv = f"{CYCLE_DIR}/seed_{seed}/esm3_outputs/esm3_scores.csv"
    if not os.path.exists(esm3_csv): continue
    seen = {}
    with open(esm3_csv) as f:
        for row in csv.DictReader(f):
            orig_id = extract_orig_id(row["id"])
            if not orig_id: continue
            try:
                ptm = float(row["ptm"])
                if orig_id not in seen or ptm > seen[orig_id]: seen[orig_id] = ptm
            except: pass
    for nid,(s,oid) in id_mapping.items():
        if s==seed and oid in seen: esm3[nid]=seen[oid]
print(f"  ESM3: {len(esm3)}")

# Chai
print("Extracting Chai...")
chai_pl = {}
for npz in glob.glob(f"{CYCLE_DIR}/chai/outputs/*/scores.model_idx_*.npz"):
    did = os.path.basename(os.path.dirname(npz))
    m = re.search(r'_id(\d+)$', did)
    if not m: continue
    nid = m.group(1)
    d = np.load(npz, allow_pickle=True)
    if "per_chain_pair_iptm" not in d: continue
    mat = d["per_chain_pair_iptm"]
    if mat.ndim == 3: mat = mat[0]
    if ci_fmn < mat.shape[0]:
        s = (mat[pi,ci_fmn]+mat[ci_fmn,pi])/2
        if nid not in chai_pl or s > chai_pl[nid]: chai_pl[nid]=float(s)
print(f"  Chai: {len(chai_pl)}")

# AF3
print("Extracting AF3...")
af3 = {}
af3_cif_paths = {}
for jf in glob.glob(f"{CYCLE_DIR}/af3/outputs/batch_*/id*/*summary_confidences.json"):
    parts = jf.split('/')
    if 'seed-' in parts[-2]: continue
    nid = re.search(r'id(\d+)$', parts[-2])
    if not nid: continue
    nid = nid.group(1)
    with open(jf) as f: d = json.load(f)
    pc = d.get("chain_pair_iptm",[])
    if len(pc) > ci_fmn: af3[nid] = float(pc[pi][ci_fmn])
    cif = glob.glob(os.path.join(os.path.dirname(jf),"*_model.cif"))
    if cif: af3_cif_paths[nid] = cif[0]
print(f"  AF3: {len(af3)} CIFs: {len(af3_cif_paths)}")

# Burial and gate
print("Computing burial and gate...")
burial = {}
coord_gate = {}
parser = MMCIFParser(QUIET=True)
for nid, cif_path in af3_cif_paths.items():
    try:
        structure = parser.get_structure("s", cif_path)[0]
        protein_atoms = np.array([a.coord for a in structure["A"].get_atoms()])
        fmn_atoms = np.array([a.coord for a in structure["C"].get_atoms() if a.get_parent().resname.strip()=="FMN"])
        hem_atoms = np.array([a.coord for a in structure["B"].get_atoms() if a.get_parent().resname.strip()=="HEM"])
        if len(fmn_atoms)==0 or len(hem_atoms)==0: continue
        fmn_burial = sum(1 for p in protein_atoms if np.linalg.norm(fmn_atoms-p,axis=1).min()<=8.0)
        hem_burial = sum(1 for p in protein_atoms if np.linalg.norm(hem_atoms-p,axis=1).min()<=8.0)
        burial[nid] = fmn_burial/hem_burial if hem_burial>0 else 0.0
        fe = his9 = his67 = None
        for res in structure["B"]:
            if res.resname.strip()=="HEM" and "FE" in res: fe=res["FE"].coord
        for res in structure["A"]:
            if res.id[1]==9 and "NE2" in res: his9=res["NE2"].coord
            if res.id[1]==67 and "NE2" in res: his67=res["NE2"].coord
        if fe is not None and his9 is not None and his67 is not None:
            coord_gate[nid] = bool(np.linalg.norm(his9-fe)<3.0 and np.linalg.norm(his67-fe)<3.0)
        else: coord_gate[nid] = False
    except: pass
print(f"  Burial: {len(burial)} Gate: {len(coord_gate)}")

# Compile
all_ids = set(chai_pl)|set(af3)
results = []
for nid in all_ids:
    c=chai_pl.get(nid,0.0); a=af3.get(nid,0.0)
    b=boltz_pl.get(nid,0.0); e=esm3.get(nid,0.0)
    bur=burial.get(nid,0.0); gate=coord_gate.get(nid,False)
    lineage,orig_id=id_mapping.get(nid,("unknown",nid))
    results.append({
        "id":nid,"lineage":lineage,"original_id":orig_id,
        "boltz_pl":round(b,4),"esm3_ptm":round(e,4),
        "chai_pl":round(c,4),"af3_pl":round(a,4),
        "burial":round(bur,4),"haem_coordination_gate":gate,
        "track7_all5_mean":round((a+c+b+e+bur)/5,4)
    })

out_csv = f"{TRAJ}/all_scores_FMN_cycle{cycle}.csv"
fn = ["id","lineage","original_id","boltz_pl","esm3_ptm","chai_pl","af3_pl","burial","haem_coordination_gate","track7_all5_mean"]
with open(out_csv,"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fn); w.writeheader(); w.writerows(results)

gated = [r for r in results if r["haem_coordination_gate"]]
print(f"\n{len(gated)} of {len(results)} pass gate")
gated.sort(key=lambda r: -r["track7_all5_mean"])
print("\nTop 10 gated:")
for r in gated[:10]:
    print(f"  id{r['id']} ({r['lineage']} orig{r['original_id']}): score={r['track7_all5_mean']:.4f} boltz={r['boltz_pl']:.3f} esm3={r['esm3_ptm']:.3f} chai={r['chai_pl']:.3f} af3={r['af3_pl']:.3f} burial={r['burial']:.3f}")
print(f"\nSaved to {out_csv}")
