# RQ2 step 7 compile for next seed selection - FMN
# This one adds the burial ratio and haem-coordination gate ontop of the normal multi track compile 
# Chain layout here: A=protein, B=retained HEM (Hem1 in this track's naming), C=FMN.
import json, csv, glob, os, re
import numpy as np
from Bio.PDB import MMCIFParser

ci = 2
pi = 0

chai_output_path = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_FMN_hemB/cycle_1/chai/outputs"
af3_output_path = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_FMN_hemB/cycle_1/af3/outputs"
boltz_predictions = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_FMN_hemB/cycle_1/boltz/outputs/boltz_results_input/predictions"
esm3_csv = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_FMN_hemB/cycle_1/esm3/outputs/esm3_scores.csv"
trajectory_dir = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_FMN_hemB/design_trajectory/cycle_1"
os.makedirs(trajectory_dir, exist_ok=True)

def extract_num(s): ## tries standard 'idN' pattern first, falls back to ligMPNN packed_N_ naming if it doesnt match
    m = re.search(r'(?:^|_)id(\d+)(?:$|_)', s)
    if m:
        return m.group(1)
    m = re.search(r'packed_(\d+)_', s)
    if m:
        return m.group(1)
    return s

# boltz
boltz_pl = {}
boltz_conf = {}
for seq_dir in glob.glob(f"{boltz_predictions}/*/"):
    sid_raw = os.path.basename(seq_dir.rstrip('/'))
    sid = extract_num(sid_raw)
    vals = []
    confs = []
    for jf in glob.glob(f"{seq_dir}confidence_*.json"):
        with open(jf) as f: d = json.load(f)
        confs.append(d.get("confidence_score", 0.0))
        pc = d.get("pair_chains_iptm")
        if pc:
            try: vals.append((pc[str(pi)][str(ci)] + pc[str(ci)][str(pi)])/2)
            except (KeyError, TypeError): pass
    if vals: boltz_pl[sid] = max(vals)
    if confs: boltz_conf[sid] = sum(confs)/len(confs) ## collected but unused 

# chai
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

# esm3
esm3 = {}
if os.path.exists(esm3_csv):
    with open(esm3_csv) as f:
        for row in csv.DictReader(f):
            sid = extract_num(row["id"])
            try: esm3[sid] = float(row["ptm"])
            except (KeyError, ValueError): pass

# af3
af3 = {}
af3_cif_paths = {}
files = glob.glob(f"{af3_output_path}/batch_*/*_summary_confidences.json") + \
        glob.glob(f"{af3_output_path}/batch_*/*/*_summary_confidences.json")
for f in files:
    try:
        with open(f) as fh: d = json.load(fh)
        sid = extract_num(os.path.basename(f).split("_summary")[0])
        af3[sid] = d["chain_pair_iptm"][pi][ci]
        cif_candidates = glob.glob(os.path.join(os.path.dirname(f), "*_model.cif"))
        if cif_candidates:
            af3_cif_paths[sid] = cif_candidates[0]
    except Exception: pass

# Track 6: burial ratio, computed directly from the AF3 structure
# Haem coordination gate: His37/His95 must both sit under 3A from the retained haem Fe, otherwise the AF3 structure has broken haem geometry and the design is excluded from ranking
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
        his37_ne2 = None
        his95_ne2 = None
        for res in structure["A"]:
            if res.id[1] == 37 and "NE2" in res:
                his37_ne2 = res["NE2"].coord
            if res.id[1] == 95 and "NE2" in res:
                his95_ne2 = res["NE2"].coord
        if fe_coord is not None and his37_ne2 is not None and his95_ne2 is not None:
            d37 = np.linalg.norm(his37_ne2 - fe_coord)
            d95 = np.linalg.norm(his95_ne2 - fe_coord)
            coord_gate[sid] = bool(d37 < 3.0 and d95 < 3.0)
        else:
            coord_gate[sid] = False
    except Exception:
        pass

all_ids = set(chai_pl) | set(af3)
print(f"Boltz ids: {len(boltz_pl)}, ESM3 ids: {len(esm3)}, Chai ids: {len(chai_pl)}, AF3 ids: {len(af3)}, Burial ids: {len(burial)}, union(chai,af3): {len(all_ids)}")

results = []
for sid in all_ids:
    c = chai_pl.get(sid, 0.0); a = af3.get(sid, 0.0)
    b = boltz_pl.get(sid, 0.0); e = esm3.get(sid, 0.0)
    bur = burial.get(sid, 0.0)
    gate = coord_gate.get(sid, False)
    results.append({"id": sid, "boltz_pl": b, "esm3_ptm": e, "chai_pl": c, "af3_pl": a, "burial_ratio": bur,
                    "haem_coordination_gate": gate,
                    "track1_af3only": a, "track2_af3chai": (a+c)/2,
                    "track3_all4": (a+c+b+e)/4, "track6_burial": bur,
                    "track7_all5_mean": (a+c+b+e+bur)/5})
ar = sorted(results, key=lambda r: -r["af3_pl"]); cr = sorted(results, key=lambda r: -r["chai_pl"])
arank = {r["id"]: i+1 for i,r in enumerate(ar)}; crank = {r["id"]: i+1 for i,r in enumerate(cr)}
for r in results:
    r["track4_rank_sum"] = arank[r["id"]] + crank[r["id"]]
    r["track5_chai_only"] = r["chai_pl"]

out_csv = f"{trajectory_dir}/all_scores_FMN_hem1_v2.csv"
fn = ["id","boltz_pl","esm3_ptm","chai_pl","af3_pl","burial_ratio","haem_coordination_gate","track1_af3only","track2_af3chai",
      "track3_all4","track4_rank_sum","track5_chai_only","track6_burial","track7_all5_mean"]
with open(out_csv,"w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=fn); w.writeheader(); w.writerows(results)
print(f"Scores for {len(results)} sequences written to {out_csv}")

specs = [("track1_af3only","Track1_AF3",True),("track2_af3chai","Track2_AF3Chai_mean",True),
         ("track3_all4","Track3_All4_mean",True),("track4_rank_sum","Track4_RankSumConsensus",False),
         ("track5_chai_only","Track5_Chai_only",True),("track6_burial","Track6_Burial",True),
         ("track7_all5_mean","Track7_All5_mean",True)]
seeds_csv = f"{trajectory_dir}/next_cycle_shortlist_v2.csv"
seeds_txt = f"{trajectory_dir}/next_cycle_shortlist_v2.txt"
gated_results = [r for r in results if r["haem_coordination_gate"]]
print(f"\n{len(gated_results)} of {len(results)} designs pass the haem coordination gate (His37 and His95 both under 3A)")

seed_rows = [] ## this script keeps top 3 shortlist per track, which are then expert-in-the-loop inspected 
lines = [f"RQ2 FMN-at-HEM_B cycle 1 -- shortlist (top 3 per track, haem-coordination-gated)",
         f"From {len(gated_results)} of {len(results)} scored sequences passing the gate", ""]
for k, nm, hib in specs:
    if not gated_results:
        lines.append(f"{nm}: NO RESULTS PASS GATE")
        continue
    ranked = sorted(gated_results, key=lambda r: r[k], reverse=hib)
    top3 = ranked[:3]
    lines.append(f"{nm}:")
    for rank, cand in enumerate(top3, start=1):
        seed_rows.append({"track": nm, "rank": rank, "id": cand["id"], "track_score": cand[k],
                          "boltz_pl": cand["boltz_pl"], "esm3_ptm": cand["esm3_ptm"],
                          "chai_pl": cand["chai_pl"], "af3_pl": cand["af3_pl"], "burial_ratio": cand["burial_ratio"],
                          "haem_coordination_gate": cand["haem_coordination_gate"]})
        lines.append(f"  {rank}. id{cand['id']} (score={cand[k]:.4f}, boltz={cand['boltz_pl']:.3f}, "
                     f"esm3={cand['esm3_ptm']:.3f}, chai={cand['chai_pl']:.3f}, af3={cand['af3_pl']:.3f}, burial={cand['burial_ratio']:.3f})")
with open(seeds_csv,"w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=["track","rank","id","track_score","boltz_pl","esm3_ptm","chai_pl","af3_pl","burial_ratio","haem_coordination_gate"])
    w.writeheader(); w.writerows(seed_rows)

uniq = sorted(set(r["id"] for r in seed_rows), key=lambda x: int(x)) if seed_rows else []
lines.append("")
lines.append(f"Unique seed ids ({len(uniq)}): " + ", ".join(f"id{x}" for x in uniq))
with open(seeds_txt,"w") as f: f.write("\n".join(lines)+"\n")
for l in lines: print(l)
