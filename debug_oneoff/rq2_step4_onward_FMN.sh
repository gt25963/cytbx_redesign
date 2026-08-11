#!/usr/bin/env bash
#Exit if any referenced variable is undefined (catches typos early)
set -u

# Step 4: score Boltz outputs on PROTEIN<->COFACTOR pair [0][2]. All variables here (biopython_python, boltz_predictions, cofactor, etc.) are exported by resume_rq2_step4.sh before this script is sourced/called
echo "step 4: scoring Boltz-2 (protein<->cofactor pair) at $(timestamp)"
"${biopython_python}" - << EOF
import json, csv, glob, os
ci = ${cofactor_chain_index}     # cofactor chain index (2)
pi = 0                            # protein chain index
predictions_dir = "${boltz_predictions}"
output_csv = "${boltz_output_path}/scores_file.csv"
rows = []
#Extract the protein<->cofactor pairwise iptm for every candidate, averaging
#both directions and taking the best model per candidate
for seq_dir in glob.glob(f"{predictions_dir}/*/"):
    seq_id = os.path.basename(seq_dir.rstrip('/'))
    json_files = glob.glob(f"{seq_dir}confidence_*.json")
    if not json_files: continue
    conf_scores, pl_scores = [], []
    for jf in json_files:
        with open(jf) as f: d = json.load(f)
        conf_scores.append(d.get("confidence_score", 0.0))
        pc = d.get("pair_chains_iptm")
        if pc:
            try:
                v = (pc[str(pi)][str(ci)] + pc[str(ci)][str(pi)]) / 2
                pl_scores.append(v)
            except (KeyError, TypeError):
                pass
    mc = sum(conf_scores)/len(conf_scores) if conf_scores else 0.0
    mp = sum(pl_scores)/len(pl_scores) if pl_scores else 0.0
    rows.append({"id": seq_id, "confidence_score": mc, "iptm": mp,
                 "calculated_average": (mc+mp)/2})
rows.sort(key=lambda x: x["calculated_average"], reverse=True)
with open(output_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["id","confidence_score","iptm","calculated_average"])
    w.writeheader(); w.writerows(rows)
print(f"scored {len(rows)} sequences (protein<->cofactor iptm)")
EOF
cp "${boltz_output_path}/scores_file.csv" "${trajectory_path}/cycle_${i}/"

#Combine Boltz-2 and ESM3 scores into one file, used to rank and pick the top 50
echo "step 4b: combine Boltz + ESM3 at $(timestamp)"
"${biopython_python}" "${work_directory}/scripts/combine_scores.py" \
    "${boltz_output_path}/scores_file.csv" \
    "${esm3_output_path}/esm3_scores.csv" \
    "${boltz_output_path}/combined_scores.csv"
cp "${boltz_output_path}/combined_scores.csv" "${trajectory_path}/cycle_${i}/"

#Read off the top N candidates (Boltz+ESM3 fast filter) to carry forward to Chai-1/AF3
top_ids=()
while IFS=',' read -r id rest; do
    [ "${id}" == "id" ] && continue
    top_ids+=("${id}")
    [ "${#top_ids[@]}" -ge "${top_n_for_chai_af3}" ] && break
done < "${boltz_output_path}/combined_scores.csv"
echo "top ${#top_ids[@]} selected for Chai-1 + AF3: ${top_ids[*]}"

numeric_ids=()
for id in "${top_ids[@]}"; do
    num=$(echo "${id}" | grep -oP '\d+$'); [ -n "${num}" ] && numeric_ids+=("${num}")
done

#  Step 5: Chai-1 array (single protein chain + 2 ligands) 
echo "step 5: Chai-1 on top ${#numeric_ids[@]} at $(timestamp)"
#Pull out just the top-50 sequences' fasta entries, then convert to Chai-1's multi-record format (protein + retained haem + cofactor)
top_ids_pattern=$(printf "| id=%s," "${numeric_ids[@]}"); top_ids_pattern="${top_ids_pattern:2}"
grep -A1 -E "(${top_ids_pattern})" "${fasta_file}" > "${chai_input_path}/top_sequences.fa" 2>/dev/null || true
"${biopython_python}" "${work_directory}/scripts/prep_chai_fasta_rq2.py" \
    "${chai_input_path}/top_sequences.fa" \
    "${chai_input_path}/chai_input.fa" \
    "${cofactor}"
chai_ids_file="${chai_input_path}/chai_ids.txt"; > "${chai_ids_file}"
for id in "${numeric_ids[@]}"; do echo "${id}" >> "${chai_ids_file}"; done
n_chai_ids=$(wc -l < "${chai_ids_file}")
sbatch --array=0-$((n_chai_ids - 1))%5 "${work_directory}/submit_chai_array.sh" \
    "${chai_ids_file}" "${chai_input_path}/chai_input.fa" "${chai_output_path}"

#Poll until the Chai-1 array finishes
while true; do
    n_chai_jobs=$(squeue -u "$(whoami)" -n chai_array 2>/dev/null | wc -l)
    [ "${n_chai_jobs}" -le 1 ] && { echo "Chai-1 no longer queued, verifying..."; break; }
    echo "Chai-1 array running... ${n_chai_jobs} remaining (incl header)"; sleep 60
done
#Verify actual output files exist, don't just trust SLURM queue status
missing_chai_ids=()
for id in "${numeric_ids[@]}"; do
    find "${chai_output_path}" -path "*id${id}*combined_scores.csv" 2>/dev/null | grep -q . || missing_chai_ids+=("${id}")
done
if [ "${#missing_chai_ids[@]}" -gt 0 ]; then
    echo "WARNING: ${#missing_chai_ids[@]} Chai-1 missing: ${missing_chai_ids[*]}"
else
    echo "Chai-1 complete: all ${n_chai_ids} present"
fi

#  Step 6: AF3 array 
echo "step 6: AF3 (array) on top ${#numeric_ids[@]} at $(timestamp)"
"${biopython_python}" "${work_directory}/scripts/fasta_to_af3_nomsa_rq2.py" \
    "${chai_input_path}/top_sequences.fa" "${af3_input_path}" "${cofactor}"
#Split into batches of 5 for the array job
cd "${af3_input_path}"; mkdir -p batches; batch_i=0; batch=0
for f in id*.json; do
    [ $((batch_i % 5)) -eq 0 ] && { batch=$((batch+1)); mkdir -p "batches/batch_${batch}"; }
    mv "${f}" "batches/batch_${batch}/"; batch_i=$((batch_i+1))
done
n_batches=$(ls "${af3_input_path}/batches" | wc -l); cd "${work_directory}"
af3_job_id=$(sbatch --parsable --array=1-${n_batches} "${work_directory}/run_af3_array.sh" \
    "${exec_directory}/cycle_${i}/af3")
echo "AF3 array job: ${af3_job_id}"
#Poll until the AF3 array finishes
while true; do
    n_done=$(squeue -j "${af3_job_id}" 2>/dev/null | wc -l)
    [ "${n_done}" -le 1 ] && { echo "AF3 no longer queued, verifying..."; break; }
    echo "AF3 array running... ${n_done} remaining (incl header)"; sleep 60
done
#Verify actual output files, not just queue status
missing_af3_ids=()
for id in "${numeric_ids[@]}"; do
    find "${af3_output_path}" -path "*id${id}*summary_confidences.json" ! -path "*seed*" 2>/dev/null | grep -q . || missing_af3_ids+=("${id}")
done
if [ "${#missing_af3_ids[@]}" -gt 0 ]; then
    echo "WARNING: ${#missing_af3_ids[@]} AF3 missing: ${missing_af3_ids[*]}"
else
    echo "AF3 complete: all ${#numeric_ids[@]} present"
fi

# Step 7: compile + five-track selection (protein<->cofactor [0][2]) 
echo "step 7: compiling scores + seed selection at $(timestamp)"
"${biopython_python}" - << EOF
import json, csv, glob, os
import numpy as np
ci = ${cofactor_chain_index}; pi = 0

#Boltz-2: re-extract true protein-cofactor pairwise score for the top-50 (same corrected extraction logic as Step 4, kept best model per candidate)
boltz_predictions = "${boltz_predictions}"
boltz_pl = {}
for seq_dir in glob.glob(f"{boltz_predictions}/*/"):
    sid = os.path.basename(seq_dir.rstrip('/'))
    vals = []
    for jf in glob.glob(f"{seq_dir}confidence_*.json"):
        with open(jf) as f: d = json.load(f)
        pc = d.get("pair_chains_iptm")
        if pc:
            try: vals.append((pc[str(pi)][str(ci)] + pc[str(ci)][str(pi)])/2)
            except (KeyError, TypeError): pass
    if vals: boltz_pl[sid] = max(vals)

#Chai-1: read the corrected protein-cofactor pairwise score directly from raw npz
chai_pl = {}
for npz in glob.glob("${chai_output_path}/*/scores.model_idx_*.npz"):
    sid = os.path.basename(os.path.dirname(npz)).replace("top_scoring.cif_","").replace("id","")
    d = np.load(npz)
    if "per_chain_pair_iptm" not in d: continue
    m = d["per_chain_pair_iptm"]
    if m.ndim == 3: m = m[0]
    if ci < m.shape[0] and pi < m.shape[1]:
        s = (m[pi, ci] + m[ci, pi]) / 2
        if sid not in chai_pl or s > chai_pl[sid]: chai_pl[sid] = s

#ESM3 pTM, from the compiled scores CSV for this cycle
esm3 = {}
ec = "${esm3_output_path}/esm3_scores.csv"
if os.path.exists(ec):
    with open(ec) as f:
        for row in csv.DictReader(f):
            sid = row["id"].replace("top_scoring.cif_","").replace("id","")
            try: esm3[sid] = float(row["ptm"])
            except (KeyError, ValueError): pass

#AF3: protein-cofactor pairwise score from chain_pair_iptm
af3 = {}
files = glob.glob("${af3_output_path}/batch_*/*_summary_confidences.json") + \
        glob.glob("${af3_output_path}/batch_*/*/*_summary_confidences.json")
for f in files:
    try:
        with open(f) as fh: d = json.load(fh)
        sid = os.path.basename(f).split("_summary")[0].replace("id","")
        af3[sid] = d["chain_pair_iptm"][pi][ci]
    except Exception: pass

#Combine all four tools' scores into one row per candidate, computing the five selection tracks used for seed advancement each cycle
all_ids = set(chai_pl) | set(af3)
results = []
for sid in all_ids:
    c = chai_pl.get(sid, 0.0); a = af3.get(sid, 0.0)
    b = boltz_pl.get(sid, 0.0); e = esm3.get(sid, 0.0)
    results.append({"id": sid, "boltz_pl": b, "esm3_ptm": e, "chai_pl": c, "af3_pl": a,
                    "track1_af3only": a, "track2_af3chai": (a+c)/2,
                    "track3_all4": (a+c+b+e)/4})
#Track4 = combined rank-sum of AF3 and Chai-1 rankings (lower is better)
ar = sorted(results, key=lambda r: -r["af3_pl"]); cr = sorted(results, key=lambda r: -r["chai_pl"])
arank = {r["id"]: i+1 for i,r in enumerate(ar)}; crank = {r["id"]: i+1 for i,r in enumerate(cr)}
for r in results:
    r["track4_rank_sum"] = arank[r["id"]] + crank[r["id"]]
    r["track5_chai_only"] = r["chai_pl"]

#Write the full per-candidate scores table for this cycle (feeds Table 3)
out_csv = "${trajectory_path}/cycle_${i}/all_scores_${cofactor}.csv"
fn = ["id","boltz_pl","esm3_ptm","chai_pl","af3_pl","track1_af3only","track2_af3chai",
      "track3_all4","track4_rank_sum","track5_chai_only"]
with open(out_csv,"w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=fn); w.writeheader(); w.writerows(results)
print(f"Scores for {len(results)} sequences -> {out_csv}")

#Report the winning candidate under each of the five selection tracks, and write a summary of unique seed candidates advancing to the next cycle
specs = [("track1_af3only","Track1_AF3",True),("track2_af3chai","Track2_AF3Chai_mean",True),
         ("track3_all4","Track3_All4_mean",True),("track4_rank_sum","Track4_RankSumConsensus",False),
         ("track5_chai_only","Track5_Chai_only",True)]
seeds_csv = "${trajectory_path}/cycle_${i}/next_cycle_seeds.csv"
seeds_txt = "${trajectory_path}/cycle_${i}/next_cycle_seeds.txt"
seed_rows = []; lines = [f"RQ2 ${cofactor} cycle ${i} -- per-track seed selection",
                         f"From {len(results)} scored sequences",""]
for k, nm, hib in specs:
    if not results: lines.append(f"{nm}: NO RESULTS"); continue
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
lines.append(""); lines.append(f"Unique seed ids ({len(uniq)}): " + ", ".join(f"id{x}" for x in uniq))
with open(seeds_txt,"w") as f: f.write("\n".join(lines)+"\n")
for l in lines: print(l)
EOF

echo "RQ2 ${cofactor} cycle ${i} complete at $(timestamp)"
echo "Scores: ${trajectory_path}/cycle_${i}/all_scores_${cofactor}.csv"
