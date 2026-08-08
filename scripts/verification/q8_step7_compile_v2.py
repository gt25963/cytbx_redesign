import json, csv, glob, os, re
import numpy as np
from Bio.PDB import MMCIFParser

pi = 0
ci_uq8 = 3

af3_output_path = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_Q8/cycle_1/af3/outputs"
chai_output_path = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_Q8/cycle_1/chai/outputs"
boltz_predictions = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_Q8/cycle_1/boltz/outputs/boltz_results_input/predictions"
esm3_csv = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_Q8/cycle_1/esm3/outputs/esm3_scores.csv"
trajectory_dir = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_Q8/design_trajectory/cycle_1"
os.makedirs(trajectory_dir, exist_ok=True)

def extract_num(s):
    m = re.search(r'(?:^|_)id(\d+)(?:$|_)', s)
    if m:
        return m.group(1)
    return s

boltz_pl = {}
for seq_dir in glob.glob(f"{boltz_predictions}/*/"):
    sid = extract_num(os.path.basename(seq_dir.rstrip('/')))
    vals = []
    for jf in glob.glob(f"{seq_dir}confidence_*.json"):
        with open(jf) as f: d = json.load(f)
        pc = d.get("pair_chains_iptm")
        if pc:
            try: vals.append((pc[str(pi)][str(ci_uq8)] + pc[str(ci_uq8)][str(pi)])/2)
            except (KeyError, TypeError): pass
    if vals: boltz_pl[sid] = max(vals)

chai_pl = {}
for npz in glob.glob(f"{chai_output_path}/*/scores.model_idx_*.npz"):
    sid = extract_num(os.path.basename(os.path.dirname(npz)))
    d = np.load(npz)
    if "per_chain_pair_iptm" not in d: continue
    m = d["per_chain_pair_iptm"]
    if m.ndim == 3: m = m[0]
    if ci_uq8 < m.shape[0] and pi < m.shape[1]:
        s = (m[pi, ci_uq8] + m[ci_uq8, pi]) / 2
        if sid not in chai_pl or s > chai_pl[sid]: chai_pl[sid] = s

esm3 = {}
if os.path.exists(esm3_csv):
    with open(esm3_csv) as f:
        for row in csv.DictReader(f):
            sid = extract_num(row["id"])
            try: esm3[sid] = float(row["ptm"])
            except (KeyError, ValueError): pass

af3 = {}
af3_cif_paths = {}
files = glob.glob(f"{af3_output_path}/batch_*/*_summary_confidences.json") + \
        glob.glob(f"{af3_output_path}/batch_*/*/*_summary_confidences.json")
for f in files:
    try:
        with open(f) as fh: d = json.load(fh)
        sid = extract_num(os.path.basename(f).split("_summary")[0])
        af3[sid] = d["chain_pair_iptm"][pi][ci_uq8]
        cif_candidates = glob.glob(os.path.join(os.path.dirname(f), "*_model.cif"))
        if cif_candidates:
            af3_cif_paths[sid] = cif_candidates[0]
    except Exception: pass

# Burial (Q8 vs each retained haem) and haem coordination gate (His37/His95 on HEM_B, His9/His67 on HEM_C)
burial_vs_hemb = {}
burial_vs_hemc = {}
coord_gate = {}
parser = MMCIFParser(QUIET=True)
for sid, cif_path in af3_cif_paths.items():
    try:
        structure = parser.get_structure("s", cif_path)[0]
        protein_atoms = np.array([a.coord for a in structure["A"].get_atoms()])
        uq8_atoms = np.array([a.coord for a in structure["D"].get_atoms() if a.get_parent().resname.strip() == "UQ8"])
        hemb_atoms = np.array([a.coord for a in structure["B"].get_atoms() if a.get_parent().resname.strip() == "HEM"])
        hemc_atoms = np.array([a.coord for a in structure["C"].get_atoms() if a.get_parent().resname.strip() == "HEM"])

        if len(uq8_atoms) == 0 or len(hemb_atoms) == 0 or len(hemc_atoms) == 0:
            continue

        uq8_burial = sum(1 for p in protein_atoms if np.linalg.norm(uq8_atoms - p, axis=1).min() <= 8.0)
        hemb_burial = sum(1 for p in protein_atoms if np.linalg.norm(hemb_atoms - p, axis=1).min() <= 8.0)
        hemc_burial = sum(1 for p in protein_atoms if np.linalg.norm(hemc_atoms - p, axis=1).min() <= 8.0)
        burial_vs_hemb[sid] = uq8_burial / hemb_burial if hemb_burial > 0 else 0.0
        burial_vs_hemc[sid] = uq8_burial / hemc_burial if hemc_burial > 0 else 0.0

        fe_b = None
        for res in structure["B"]:
            if res.resname.strip() == "HEM" and "FE" in res:
                fe_b = res["FE"].coord
                break
        fe_c = None
        for res in structure["C"]:
            if res.resname.strip() == "HEM" and "FE" in res:
                fe_c = res["FE"].coord
                break

        his37_ne2 = his95_ne2 = his9_ne2 = his67_ne2 = None
        for res in structure["A"]:
            if res.id[1] == 37 and "NE2" in res: his37_ne2 = res["NE2"].coord
            if res.id[1] == 95 and "NE2" in res: his95_ne2 = res["NE2"].coord
            if res.id[1] == 9 and "NE2" in res: his9_ne2 = res["NE2"].coord
            if res.id[1] == 67 and "NE2" in res: his67_ne2 = res["NE2"].coord

        gate_b = False
        gate_c = False
        if fe_b is not None and his37_ne2 is not None and his95_ne2 is not None:
            gate_b = bool(np.linalg.norm(his37_ne2 - fe_b) < 3.0 and np.linalg.norm(his95_ne2 - fe_b) < 3.0)
        if fe_c is not None and his9_ne2 is not None and his67_ne2 is not None:
            gate_c = bool(np.linalg.norm(his9_ne2 - fe_c) < 3.0 and np.linalg.norm(his67_ne2 - fe_c) < 3.0)

        coord_gate[sid] = bool(gate_b and gate_c)
    except Exception:
        pass

all_ids = set(chai_pl) | set(af3)
print(f"Boltz ids: {len(boltz_pl)}, ESM3 ids: {len(esm3)}, Chai ids: {len(chai_pl)}, AF3 ids: {len(af3)}, Burial ids: {len(burial_vs_hemc)}, union(chai,af3): {len(all_ids)}")

results = []
for sid in all_ids:
    c = chai_pl.get(sid, 0.0); a = af3.get(sid, 0.0)
    b = boltz_pl.get(sid, 0.0); e = esm3.get(sid, 0.0)
    bur_hemb = burial_vs_hemb.get(sid, 0.0)
    bur_hemc = burial_vs_hemc.get(sid, 0.0)
    gate = coord_gate.get(sid, False)
    results.append({"id": sid, "boltz_pl": b, "esm3_ptm": e, "chai_pl": c, "af3_pl": a,
                    "burial_vs_hemb": bur_hemb, "burial_vs_hemc": bur_hemc,
                    "haem_coordination_gate": gate,
                    "track1_af3only": a, "track2_af3chai": (a+c)/2,
                    "track3_all4": (a+c+b+e)/4, "track6_burial": bur_hemc,
                    "track7_all5_mean": (a+c+b+e+bur_hemc)/5})

ar = sorted(results, key=lambda r: -r["af3_pl"]); cr = sorted(results, key=lambda r: -r["chai_pl"])
arank = {r["id"]: i+1 for i,r in enumerate(ar)}; crank = {r["id"]: i+1 for i,r in enumerate(cr)}
for r in results:
    r["track4_rank_sum"] = arank[r["id"]] + crank[r["id"]]
    r["track5_chai_only"] = r["chai_pl"]

out_csv = f"{trajectory_dir}/all_scores_Q8_v2.csv"
fn = ["id","boltz_pl","esm3_ptm","chai_pl","af3_pl","burial_vs_hemb","burial_vs_hemc","haem_coordination_gate",
      "track1_af3only","track2_af3chai","track3_all4","track4_rank_sum","track5_chai_only","track6_burial","track7_all5_mean"]
with open(out_csv,"w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=fn); w.writeheader(); w.writerows(results)
print(f"Scores for {len(results)} sequences written to {out_csv}")

gated_results = [r for r in results if r["haem_coordination_gate"]]
print(f"\n{len(gated_results)} of {len(results)} designs pass the haem coordination gate (His37/His95 to HEM_B Fe AND His9/His67 to HEM_C Fe, all under 3A)")

specs = [("track1_af3only","Track1_AF3",True),("track2_af3chai","Track2_AF3Chai_mean",True),
         ("track3_all4","Track3_All4_mean",True),("track4_rank_sum","Track4_RankSumConsensus",False),
         ("track5_chai_only","Track5_Chai_only",True),("track6_burial","Track6_Burial_vs_HEMC",True),
         ("track7_all5_mean","Track7_All5_mean",True)]

seeds_csv = f"{trajectory_dir}/next_cycle_shortlist_v2.csv"
seeds_txt = f"{trajectory_dir}/next_cycle_shortlist_v2.txt"
seed_rows = []
lines = [f"RQ2 Q8 cycle 1 -- shortlist (top 3 per track, haem-coordination-gated)",
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
                          "chai_pl": cand["chai_pl"], "af3_pl": cand["af3_pl"],
                          "burial_vs_hemc": cand["burial_vs_hemc"], "haem_coordination_gate": cand["haem_coordination_gate"]})
        lines.append(f"  {rank}. id{cand['id']} (score={cand[k]:.4f}, boltz={cand['boltz_pl']:.3f}, "
                     f"esm3={cand['esm3_ptm']:.3f}, chai={cand['chai_pl']:.3f}, af3={cand['af3_pl']:.3f}, burial_vs_hemc={cand['burial_vs_hemc']:.3f})")
with open(seeds_csv,"w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=["track","rank","id","track_score","boltz_pl","esm3_ptm","chai_pl","af3_pl","burial_vs_hemc","haem_coordination_gate"])
    w.writeheader(); w.writerows(seed_rows)

uniq = sorted(set(r["id"] for r in seed_rows), key=lambda x: int(x)) if seed_rows else []
lines.append("")
lines.append(f"Unique shortlist ids ({len(uniq)}): " + ", ".join(f"id{x}" for x in uniq))
with open(seeds_txt,"w") as f: f.write("\n".join(lines)+"\n")
for l in lines: print(l)
