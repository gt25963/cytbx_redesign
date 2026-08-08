import json, csv, glob, os, re
import numpy as np
from Bio.PDB import MMCIFParser

pi = 0
ci_fmn = 2

cycle2_dir = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_FMN_hemB/cycle_2"
af3_output_path = f"{cycle2_dir}/af3/outputs"
chai_output_path = f"{cycle2_dir}/chai_outputs"
trajectory_dir = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_FMN_hemB/design_trajectory/cycle_2"
os.makedirs(trajectory_dir, exist_ok=True)

id_mapping = {}
with open(f"{cycle2_dir}/pooled_id_mapping.csv") as f:
    for row in csv.DictReader(f):
        id_mapping[row["numeric_id"]] = (row["lineage"], row["original_id"])

def extract_num(s):
    m = re.search(r'(?:^|_)id(\d+)(?:$|_)', s)
    if m:
        return m.group(1)
    return s

boltz_pl = {}
for numeric_id, (lineage, orig_id) in id_mapping.items():
    boltz_dir = f"{cycle2_dir}/seed_{lineage}/boltz_outputs/boltz_results_input/predictions"
    matches = glob.glob(f"{boltz_dir}/*id{orig_id}/") + glob.glob(f"{boltz_dir}/*id{orig_id}_[0-9]*/")
    exact_matches = [m for m in matches if extract_num(os.path.basename(m.rstrip('/'))) == orig_id]
    if not exact_matches:
        continue
    vals = []
    for seq_dir in exact_matches:
        for jf in glob.glob(f"{seq_dir}confidence_*.json"):
            with open(jf) as f: d = json.load(f)
            pc = d.get("pair_chains_iptm")
            if pc:
                try: vals.append((pc[str(pi)][str(ci_fmn)] + pc[str(ci_fmn)][str(pi)])/2)
                except (KeyError, TypeError): pass
    if vals:
        boltz_pl[numeric_id] = max(vals)

print(f"Boltz matched: {len(boltz_pl)} of {len(id_mapping)}")

chai_pl = {}
for npz in glob.glob(f"{chai_output_path}/*/scores.model_idx_*.npz"):
    sid = extract_num(os.path.basename(os.path.dirname(npz)))
    d = np.load(npz)
    if "per_chain_pair_iptm" not in d: continue
    m = d["per_chain_pair_iptm"]
    if m.ndim == 3: m = m[0]
    if ci_fmn < m.shape[0] and pi < m.shape[1]:
        s = (m[pi, ci_fmn] + m[ci_fmn, pi]) / 2
        if sid not in chai_pl or s > chai_pl[sid]: chai_pl[sid] = s

print(f"Chai matched: {len(chai_pl)}")

esm3 = {}
for lineage in ["id111", "id93"]:
    esm3_csv = f"{cycle2_dir}/seed_{lineage}/esm3_outputs/esm3_scores.csv"
    if os.path.exists(esm3_csv):
        with open(esm3_csv) as f:
            for row in csv.DictReader(f):
                orig_id = extract_num(row["id"])
                for numeric_id, (lin, oid) in id_mapping.items():
                    if lin == lineage and oid == orig_id:
                        try: esm3[numeric_id] = float(row["ptm"])
                        except (KeyError, ValueError): pass

print(f"ESM3 matched: {len(esm3)}")

af3 = {}
af3_cif_paths = {}
files = glob.glob(f"{af3_output_path}/batch_*/*_summary_confidences.json") + \
        glob.glob(f"{af3_output_path}/batch_*/*/*_summary_confidences.json")
for f in files:
    try:
        with open(f) as fh: d = json.load(fh)
        sid = extract_num(os.path.basename(f).split("_summary")[0])
        af3[sid] = d["chain_pair_iptm"][pi][ci_fmn]
        cif_candidates = glob.glob(os.path.join(os.path.dirname(f), "*_model.cif"))
        if cif_candidates:
            af3_cif_paths[sid] = cif_candidates[0]
    except Exception: pass

print(f"AF3 matched: {len(af3)}")

burial = {}
coord_gate = {}
parser = MMCIFParser(QUIET=True)
for sid, cif_path in af3_cif_paths.items():
    try:
        structure = parser.get_structure("s", cif_path)[0]
        protein_atoms = np.array([a.coord for a in structure["A"].get_atoms()])
        fmn_atoms = np.array([a.coord for a in structure["C"].get_atoms() if a.get_parent().resname.strip() == "FMN"])
        hem_atoms = np.array([a.coord for a in structure["B"].get_atoms() if a.get_parent().resname.strip() == "HEM"])
        if len(fmn_atoms) == 0 or len(hem_atoms) == 0:
            continue
        fmn_burial = sum(1 for p in protein_atoms if np.linalg.norm(fmn_atoms - p, axis=1).min() <= 8.0)
        hem_burial = sum(1 for p in protein_atoms if np.linalg.norm(hem_atoms - p, axis=1).min() <= 8.0)
        burial[sid] = fmn_burial / hem_burial if hem_burial > 0 else 0.0

        fe_coord = None
        for res in structure["B"]:
            if res.resname.strip() == "HEM" and "FE" in res:
                fe_coord = res["FE"].coord
                break
        his37_ne2 = his95_ne2 = None
        for res in structure["A"]:
            if res.id[1] == 37 and "NE2" in res: his37_ne2 = res["NE2"].coord
            if res.id[1] == 95 and "NE2" in res: his95_ne2 = res["NE2"].coord
        if fe_coord is not None and his37_ne2 is not None and his95_ne2 is not None:
            d37 = np.linalg.norm(his37_ne2 - fe_coord)
            d95 = np.linalg.norm(his95_ne2 - fe_coord)
            coord_gate[sid] = bool(d37 < 3.0 and d95 < 3.0)
        else:
            coord_gate[sid] = False
    except Exception:
        pass

print(f"Coordination gate computed for {len(coord_gate)} designs")

all_ids = set(chai_pl) | set(af3)
results = []
for sid in all_ids:
    c = chai_pl.get(sid, 0.0); a = af3.get(sid, 0.0)
    b = boltz_pl.get(sid, 0.0); e = esm3.get(sid, 0.0)
    bur = burial.get(sid, 0.0)
    gate = coord_gate.get(sid, False)
    lineage, orig_id = id_mapping.get(sid, ("unknown", sid))
    results.append({"id": sid, "lineage": lineage, "original_id": orig_id,
                    "boltz_pl": b, "esm3_ptm": e, "chai_pl": c, "af3_pl": a,
                    "burial_ratio": bur, "haem_coordination_gate": gate,
                    "track7_all5_mean": (a+c+b+e+bur)/5})

out_csv = f"{trajectory_dir}/all_scores_FMN_hemB_cycle2.csv"
fn = ["id","lineage","original_id","boltz_pl","esm3_ptm","chai_pl","af3_pl","burial_ratio","haem_coordination_gate","track7_all5_mean"]
with open(out_csv,"w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=fn); w.writeheader(); w.writerows(results)

gated_results = [r for r in results if r["haem_coordination_gate"]]
print(f"\n{len(gated_results)} of {len(results)} designs pass the haem coordination gate")

gated_results.sort(key=lambda r: -r["track7_all5_mean"])
print("\nTop 5 gated, ranked by all-five-mean (incl burial):")
for r in gated_results[:5]:
    print(f"  id{r['id']} (from {r['lineage']} orig_id{r['original_id']}): score={r['track7_all5_mean']:.4f}, "
          f"boltz={r['boltz_pl']:.3f}, esm3={r['esm3_ptm']:.3f}, chai={r['chai_pl']:.3f}, af3={r['af3_pl']:.3f}, burial={r['burial_ratio']:.3f}")
