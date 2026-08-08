#!/usr/bin/env bash
# ============================================================
# CytbX RQ2 Cofactor-Swap Scoring Pipeline -- FMN variant
# Mirrors master_cytbx_4tool_C3.sh, adapted for RQ2:
#   - starts at Step 3 (designs already exist from clash-checked LigandMPNN)
#   - protein_chain_count=1; scores the PROTEIN<->COFACTOR pair [0][2]
#     (not protein-protein). Chain order: A=protein(0), B=HEM(1), C=FMN(2)
#   - retained native haem = HEM (chain B, index 1); swapped cofactor =
#     FMN (chain C, index 2). Metric of interest is pocket binding => [0][2].
#   - five-track selection retained, structurally identical to C3
# ============================================================
set -u

home_directory="/home/b5ae/mvg2713124.b5ae"
scratch_directory="/scratch/b5ae/mvg2713124.b5ae"
work_directory="${scratch_directory}/cytbx_pipeline"
master_folder="rq2/master_pipeline/RQ2_FMN"
exec_directory="${work_directory}/${master_folder}"
trajectory_path="${exec_directory}/design_trajectory"
biopython_python="${home_directory}/miniconda3/envs/biopython/bin/python"

timestamp() { date +"%F_%T"; }

# ---- RQ2 INPUTS (edit these to point at a different design set) ----
cofactor="FMN"
design_source="${work_directory}/rq2/design/FMN_pocket/cycle1_relaxed"
input_fasta="rq2/design/FMN_pocket/cycle2_combined_input.fa"
clash_csv="${design_source}/clash_check_results.csv"   # used to restrict to clash-free designs
# --------------------------------------------------------------------

ligand_spec="HEM:1,${cofactor}:1"   # retained haem + swapped cofactor, both by CCD
boltz_samples=2
top_n_for_chai_af3=50
protein_chain_count=1
cofactor_chain_index=2              # A=0 protein, B=1 HEM, C=2 cofactor
i=2

echo "RQ2 ${cofactor} pipeline starting at $(timestamp)"
for sub in boltz/inputs boltz/outputs esm3/inputs esm3/outputs \
           chai/inputs chai/outputs af3/inputs af3/outputs; do
    mkdir -p "${exec_directory}/cycle_${i}/${sub}"
done
mkdir -p "${trajectory_path}/cycle_${i}" "${work_directory}/logs"

boltz_input_path="${exec_directory}/cycle_${i}/boltz/inputs"
boltz_output_path="${exec_directory}/cycle_${i}/boltz/outputs"
boltz_predictions="${boltz_output_path}/boltz_results_yamls/predictions"
esm3_input_path="${exec_directory}/cycle_${i}/esm3/inputs"
esm3_output_path="${exec_directory}/cycle_${i}/esm3/outputs"
chai_input_path="${exec_directory}/cycle_${i}/chai/inputs"
chai_output_path="${exec_directory}/cycle_${i}/chai/outputs"
af3_input_path="${exec_directory}/cycle_${i}/af3/inputs"
af3_output_path="${exec_directory}/cycle_${i}/af3/outputs"

cp "${work_directory}/scripts/mpnn_to_boltz2_FIXED.sh" "${boltz_input_path}/mpnn_to_boltz2.sh"

# ---- Step 2.5: restrict to clash-free designs ----
echo "step 2.5: building clash-free FASTA at $(timestamp)"
clashfree_fasta="${boltz_input_path}/clashfree_input.fa"
"${biopython_python}" - << EOF
import csv, re, os
clash_csv = "${clash_csv}"
src_fa    = "${input_fasta}"
out_fa    = "${clashfree_fasta}"

# collect clash-free design numbers from the clash CSV
keep = set()
with open(clash_csv) as f:
    for row in csv.DictReader(f):
        if str(row.get("overall_clash_free")).strip().lower() == "true":
            m = re.search(r'packed_(\d+)_', row["pdb_file"])
            if m:
                keep.add(int(m.group(1)))
print(f"{len(keep)} clash-free designs to score")

# stream the FASTA, emit only records whose design id is clash-free.
# LigandMPNN headers look like: >...,  id=N, ...   (N is the design index)
kept = 0
with open(src_fa) as fin, open(out_fa, "w") as fout:
    write = False
    for line in fin:
        if line.startswith(">"):
            m = re.search(r'id=(\d+)', line)
            did = int(m.group(1)) if m else None
            write = (did in keep)
            if write:
                fout.write(line); kept += 1
        elif write:
            fout.write(line)
print(f"wrote {kept} clash-free sequences to {out_fa}")
EOF
cp "${clashfree_fasta}" "${trajectory_path}/cycle_${i}/"
fasta_file="${clashfree_fasta}"
total_sequences=$(grep -c "^>" "${fasta_file}")
echo "total clash-free sequences: ${total_sequences}"

# ---- Step 3: FASTA -> Boltz YAML (two ligands: HEM + cofactor) ----
echo "step 3: FASTA -> Boltz YAML (${ligand_spec}) at $(timestamp)"
bash "${boltz_input_path}/mpnn_to_boltz2.sh" \
    "${fasta_file}" \
    -l "${ligand_spec}" \
    -o "${boltz_input_path}/yamls"

echo "step 3b: Boltz-2 on all clash-free sequences at $(timestamp)"
sbatch "${work_directory}/submit_boltz.sh" \
    "rq2_${cofactor}_${i}" \
    "${boltz_input_path}/yamls" \
    "${boltz_output_path}" \
    "${boltz_samples}"

echo "step 3c: ESM3 on all sequences at $(timestamp)"
cp "${fasta_file}" "${esm3_input_path}/"
sbatch "${work_directory}/submit_esm3.sh" \
    "rq2_${cofactor}_${i}" \
    "${esm3_input_path}/$(basename ${fasta_file})" \
    "${esm3_output_path}"

n_yamls=$(ls "${boltz_input_path}/yamls/"*.yaml 2>/dev/null | wc -l)
total_boltz_cifs=$(( n_yamls * boltz_samples ))
while true; do
    current_cifs=$(find "${boltz_output_path}" -name "*.cif" 2>/dev/null | wc -l)
    if [ "${current_cifs}" -ge "${total_boltz_cifs}" ]; then
        echo "Boltz-2 complete at $(timestamp)"; break
    fi
    echo "Boltz-2 running... ${current_cifs}/${total_boltz_cifs}"; sleep 30
done

while true; do
    if [ -f "${esm3_output_path}/esm3_scores.csv" ]; then
        esm3_lines=$(wc -l < "${esm3_output_path}/esm3_scores.csv")
        [ "${esm3_lines}" -ge "${total_sequences}" ] && { echo "ESM3 complete at $(timestamp)"; break; }
    fi
    echo "ESM3 running..."; sleep 30
done

# ---- Step 4: score Boltz outputs on PROTEIN<->COFACTOR pair [0][2] ----
echo "step 4: scoring Boltz-2 (protein<->cofactor pair) at $(timestamp)"
"${biopython_python}" - << EOF
import json, csv, glob, os
ci = ${cofactor_chain_index}     # cofactor chain index (2)
pi = 0                            # protein chain index
predictions_dir = "${boltz_predictions}"
output_csv = "${boltz_output_path}/scores_file.csv"
rows = []
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

echo "step 4b: combine Boltz + ESM3 at $(timestamp)"
"${biopython_python}" "${work_directory}/scripts/combine_scores.py" \
    "${boltz_output_path}/scores_file.csv" \
    "${esm3_output_path}/esm3_scores.csv" \
    "${boltz_output_path}/combined_scores.csv"
cp "${boltz_output_path}/combined_scores.csv" "${trajectory_path}/cycle_${i}/"

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

# ---- Step 5: Chai-1 array (single protein chain + 2 ligands) ----
echo "step 5: Chai-1 on top ${#numeric_ids[@]} at $(timestamp)"
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

while true; do
    n_chai_jobs=$(squeue -u "$(whoami)" -n chai_array 2>/dev/null | wc -l)
    [ "${n_chai_jobs}" -le 1 ] && { echo "Chai-1 no longer queued, verifying..."; break; }
    echo "Chai-1 array running... ${n_chai_jobs} remaining (incl header)"; sleep 60
done
missing_chai_ids=()
for id in "${numeric_ids[@]}"; do
    find "${chai_output_path}" -path "*id${id}*combined_scores.csv" 2>/dev/null | grep -q . || missing_chai_ids+=("${id}")
done
if [ "${#missing_chai_ids[@]}" -gt 0 ]; then
    echo "WARNING: ${#missing_chai_ids[@]} Chai-1 missing: ${missing_chai_ids[*]}"
else
    echo "Chai-1 complete: all ${n_chai_ids} present"
fi

# ---- Step 6: AF3 array ----
echo "step 6: AF3 (array) on top ${#numeric_ids[@]} at $(timestamp)"
"${biopython_python}" "${work_directory}/scripts/fasta_to_af3_nomsa_rq2.py" \
    "${chai_input_path}/top_sequences.fa" "${af3_input_path}" "${cofactor}"
cd "${af3_input_path}"; mkdir -p batches; batch_i=0; batch=0
for f in id*.json; do
    [ $((batch_i % 5)) -eq 0 ] && { batch=$((batch+1)); mkdir -p "batches/batch_${batch}"; }
    mv "${f}" "batches/batch_${batch}/"; batch_i=$((batch_i+1))
done
n_batches=$(ls "${af3_input_path}/batches" | wc -l); cd "${work_directory}"
af3_job_id=$(sbatch --parsable --array=1-${n_batches} "${work_directory}/run_af3_array.sh" \
    "${exec_directory}/cycle_${i}/af3")
echo "AF3 array job: ${af3_job_id}"
while true; do
    n_done=$(squeue -j "${af3_job_id}" 2>/dev/null | wc -l)
    [ "${n_done}" -le 1 ] && { echo "AF3 no longer queued, verifying..."; break; }
    echo "AF3 array running... ${n_done} remaining (incl header)"; sleep 60
done
missing_af3_ids=()
for id in "${numeric_ids[@]}"; do
    find "${af3_output_path}" -path "*id${id}*summary_confidences.json" ! -path "*seed*" 2>/dev/null | grep -q . || missing_af3_ids+=("${id}")
done
if [ "${#missing_af3_ids[@]}" -gt 0 ]; then
    echo "WARNING: ${#missing_af3_ids[@]} AF3 missing: ${missing_af3_ids[*]}"
else
    echo "AF3 complete: all ${#numeric_ids[@]} present"
fi

# ---- Step 7: compile + five-track selection (protein<->cofactor [0][2]) ----
echo "step 7: compiling scores + seed selection at $(timestamp)"
"${biopython_python}" - << EOF
import json, csv, glob, os
import numpy as np
ci = ${cofactor_chain_index}; pi = 0

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

esm3 = {}
ec = "${esm3_output_path}/esm3_scores.csv"
if os.path.exists(ec):
    with open(ec) as f:
        for row in csv.DictReader(f):
            sid = row["id"].replace("top_scoring.cif_","").replace("id","")
            try: esm3[sid] = float(row["ptm"])
            except (KeyError, ValueError): pass

af3 = {}
files = glob.glob("${af3_output_path}/batch_*/*_summary_confidences.json") + \
        glob.glob("${af3_output_path}/batch_*/*/*_summary_confidences.json")
for f in files:
    try:
        with open(f) as fh: d = json.load(fh)
        sid = os.path.basename(f).split("_summary")[0].replace("id","")
        af3[sid] = d["chain_pair_iptm"][pi][ci]
    except Exception: pass

all_ids = set(chai_pl) | set(af3)
results = []
for sid in all_ids:
    c = chai_pl.get(sid, 0.0); a = af3.get(sid, 0.0)
    b = boltz_pl.get(sid, 0.0); e = esm3.get(sid, 0.0)
    results.append({"id": sid, "boltz_pl": b, "esm3_ptm": e, "chai_pl": c, "af3_pl": a,
                    "track1_af3only": a, "track2_af3chai": (a+c)/2,
                    "track3_all4": (a+c+b+e)/4})
ar = sorted(results, key=lambda r: -r["af3_pl"]); cr = sorted(results, key=lambda r: -r["chai_pl"])
arank = {r["id"]: i+1 for i,r in enumerate(ar)}; crank = {r["id"]: i+1 for i,r in enumerate(cr)}
for r in results:
    r["track4_rank_sum"] = arank[r["id"]] + crank[r["id"]]
    r["track5_chai_only"] = r["chai_pl"]

out_csv = "${trajectory_path}/cycle_${i}/all_scores_${cofactor}.csv"
fn = ["id","boltz_pl","esm3_ptm","chai_pl","af3_pl","track1_af3only","track2_af3chai",
      "track3_all4","track4_rank_sum","track5_chai_only"]
with open(out_csv,"w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=fn); w.writeheader(); w.writerows(results)
print(f"Scores for {len(results)} sequences -> {out_csv}")

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
