"""
RQ1 C3 pooled compile script.
Usage: python compile_c3_pooled.py <cycle_dir_seed1> <cycle_dir_seed2> <output_traj_dir>
Example: python compile_c3_pooled.py \
    main_pipeline/CytbX_4tool_C3_id40_cycle5/cycle_5_seed_id121 \
    main_pipeline/CytbX_4tool_C3_id40_cycle5/cycle_5_seed_id141 \
    main_pipeline/CytbX_4tool_C3_id40_cycle5/pooled_traj
"""
import json, csv, glob, os, re, sys
import numpy as np

WORK = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline"
seed_dirs = [os.path.join(WORK, sys.argv[1]), os.path.join(WORK, sys.argv[2])]
TRAJ = os.path.join(WORK, sys.argv[3])
os.makedirs(TRAJ, exist_ok=True)

protein_chain_count = 3
pairs = [(i,j) for i in range(protein_chain_count) for j in range(i+1, protein_chain_count)]

def extract_id(s):
    m = re.search(r'_id(\d+)$', s)
    return m.group(1) if m else None

all_results = []

for SEED_DIR in seed_dirs:
    seed_name = os.path.basename(SEED_DIR)

    # Boltz
    boltz = {}
    for seq_dir in glob.glob(f"{SEED_DIR}/boltz/outputs/boltz_results_input/predictions/*/"):
        sid = extract_id(os.path.basename(seq_dir.rstrip('/')))
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

    # ESM3
    esm3 = {}
    esm3_csv = f"{SEED_DIR}/esm3/outputs/esm3_scores.csv"
    if os.path.exists(esm3_csv):
        with open(esm3_csv) as f:
            for row in csv.DictReader(f):
                sid = extract_id(row["id"])
                if sid:
                    try: esm3[sid] = float(row["ptm"])
                    except: pass

    # Chai
    chai = {}
    for npz in glob.glob(f"{SEED_DIR}/chai/outputs/*/scores.model_idx_*.npz"):
        did = os.path.basename(os.path.dirname(npz))
        sid = extract_id(did)
        if not sid: continue
        d = np.load(npz, allow_pickle=True)
        if "per_chain_pair_iptm" not in d: continue
        mat = d["per_chain_pair_iptm"]
        if mat.ndim == 3: mat = mat[0]
        vals = [(mat[i][j]+mat[j][i])/2 for i,j in pairs if i < mat.shape[0] and j < mat.shape[1]]
        if vals:
            s = float(np.mean(vals))
            if sid not in chai or s > chai[sid]: chai[sid] = s

    # AF3
    af3 = {}
    for jf in glob.glob(f"{SEED_DIR}/af3/outputs/batch_*/id*/*summary_confidences.json"):
        parts = jf.split('/')
        if 'seed-' in parts[-2]: continue
        sid = re.search(r'id(\d+)$', parts[-2])
        if not sid: continue
        sid = sid.group(1)
        with open(jf) as f: d = json.load(f)
        pc = d.get("chain_pair_iptm", [])
        if len(pc) >= 3:
            vals = [(pc[i][j]+pc[j][i])/2 for i,j in pairs if i < len(pc) and j < len(pc[i])]
            if vals: af3[sid] = float(np.mean(vals))

    print(f"{seed_name}: Boltz={len(boltz)} ESM3={len(esm3)} Chai={len(chai)} AF3={len(af3)}")

    for sid in set(boltz)|set(esm3)|set(chai)|set(af3):
        b = boltz.get(sid,0.0); e = esm3.get(sid,0.0)
        c = chai.get(sid,0.0); a = af3.get(sid,0.0)
        all_results.append({
            "id": f"{seed_name}_id{sid}", "seed": seed_name, "orig_id": sid,
            "boltz_pp": round(b,4), "esm3_ptm": round(e,4),
            "chai_pp": round(c,4), "af3_pp": round(a,4),
            "combined_all4": round((b+e+c+a)/4, 4)
        })

all_results.sort(key=lambda r: -r["combined_all4"])

# Five tracks
af3_r = {r["id"]: i for i,r in enumerate(sorted(all_results, key=lambda r: -r["af3_pp"]))}
chai_r = {r["id"]: i for i,r in enumerate(sorted(all_results, key=lambda r: -r["chai_pp"]))}
for r in all_results: r["ranksum"] = af3_r[r["id"]] + chai_r[r["id"]]

tracks = {
    "Track1_AF3":     sorted(all_results, key=lambda r: -r["af3_pp"])[0],
    "Track2_AF3Chai": sorted(all_results, key=lambda r: -(r["af3_pp"]+r["chai_pp"])/2)[0],
    "Track3_All4":    sorted(all_results, key=lambda r: -r["combined_all4"])[0],
    "Track4_RankSum": sorted(all_results, key=lambda r: r["ranksum"])[0],
    "Track5_Chai":    sorted(all_results, key=lambda r: -r["chai_pp"])[0],
}

print(f"\nTop 10 by combined all-4:")
for r in all_results[:10]:
    print(f"  {r['id']}: combined={r['combined_all4']:.4f} chai={r['chai_pp']:.3f} af3={r['af3_pp']:.3f}")

print(f"\n=== Track Seeds ===")
seen = {}
for track, r in tracks.items():
    key = r["id"]
    dup = " (DUPLICATE)" if key in seen else ""
    print(f"{track}: {key} combined={r['combined_all4']:.4f} chai={r['chai_pp']:.3f}{dup}")
    if not dup: seen[key] = track

for r in all_results:
    if len(seen) >= 5: break
    if r["id"] not in seen:
        seen[r["id"]] = "Substitute"
        print(f"Substitute: {r['id']} combined={r['combined_all4']:.4f}")

print(f"\nFinal seeds: {list(seen.keys())}")

fn = ["id","seed","orig_id","boltz_pp","esm3_ptm","chai_pp","af3_pp","combined_all4","ranksum"]
with open(f"{TRAJ}/all_scores_pooled.csv","w",newline="") as f:
    w = csv.DictWriter(f,fieldnames=fn); w.writeheader(); w.writerows(all_results)
print(f"Saved to {TRAJ}/all_scores_pooled.csv")
