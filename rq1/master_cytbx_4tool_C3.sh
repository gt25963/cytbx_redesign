#!/usr/bin/env bash
# dos2unix master_cytbx_4tool_C3.sh if edited on Windows

# ============================================================
# CytbX 4-Tool Combined Pipeline -- RQ1 -- C3 Cycle 1
# Mirrors master_cytbx_4tool.sh (C2) exactly, with:
#   - starting structure: C3/holo/top3_holo.pdb (trimer)
#   - protein_chain_count=3 throughout scoring (Step 7)
#   - uses the CHAIN-COLLISION-FIXED mpnn_to_boltz2.sh
#   - Step 7 uses protein-protein-specific scores (pair_chains_iptm
#     for Boltz-2, per_chain_pair_iptm from .npz for Chai-1), not
#     the inflated all-chain-pair summaries
#   - writes next_cycle_seeds.csv / .txt as a durable per-track
#     selection artifact
#
# Residue selection for the trimer (already run today):
#   FIXED=A9 A37 A67 A95 B9 B37 B67 B95 C9 C37 C67 C95
#   INTERFACE=A66 A73 A76 A77 A79 A80 A81 A83 A84 A87 A88 A89 A90
#     A91 A92 A93 A94 A96 A97 A98 A99 A100 A101 A103 A104 A105 A108
#     A111 A113 [B.../C... mirrored]
# ============================================================

set -u

home_directory="/home/b5ae/mvg2713124.b5ae"
scratch_directory="/scratch/b5ae/mvg2713124.b5ae"
work_directory="${scratch_directory}/cytbx_pipeline"
master_folder="CytbX_4tool_C3"
exec_directory="${work_directory}/main_pipeline/${master_folder}"
trajectory_path="${exec_directory}/design_trajectory"
biopython_python="${home_directory}/miniconda3/envs/biopython/bin/python"

timestamp() { date +"%F_%T"; }

initial_structure="main_pipeline/CytbX_4tool_C3/cycle_2_seed_id35.pdb"

ligand=HEM
num_ligands=4
mpnn_model_type="global_label_membrane_mpnn"
mpnn_temp=0.5
total_sequences=200
boltz_samples=2
top_n_for_chai_af3=50
protein_chain_count=3
i=3

echo "4-tool pipeline (C3) starting at $(timestamp)"
mkdir -p "${exec_directory}/cycle_${i}/LigandMPNN/inputs"
mkdir -p "${exec_directory}/cycle_${i}/LigandMPNN/outputs"
mkdir -p "${exec_directory}/cycle_${i}/boltz/inputs"
mkdir -p "${exec_directory}/cycle_${i}/boltz/outputs"
mkdir -p "${exec_directory}/cycle_${i}/esm3/inputs"
mkdir -p "${exec_directory}/cycle_${i}/esm3/outputs"
mkdir -p "${exec_directory}/cycle_${i}/chai/inputs"
mkdir -p "${exec_directory}/cycle_${i}/chai/outputs"
mkdir -p "${exec_directory}/cycle_${i}/af3/inputs"
mkdir -p "${exec_directory}/cycle_${i}/af3/outputs"
mkdir -p "${trajectory_path}/cycle_${i}"
mkdir -p "${work_directory}/logs"

mpnn_input_path="${exec_directory}/cycle_${i}/LigandMPNN/inputs"
mpnn_output_path="${exec_directory}/cycle_${i}/LigandMPNN/outputs/seqs"
boltz_input_path="${exec_directory}/cycle_${i}/boltz/inputs"
boltz_output_path="${exec_directory}/cycle_${i}/boltz/outputs"
boltz_predictions="${boltz_output_path}/boltz_results_input/predictions"
esm3_input_path="${exec_directory}/cycle_${i}/esm3/inputs"
esm3_output_path="${exec_directory}/cycle_${i}/esm3/outputs"
chai_input_path="${exec_directory}/cycle_${i}/chai/inputs"
chai_output_path="${exec_directory}/cycle_${i}/chai/outputs"
af3_input_path="${exec_directory}/cycle_${i}/af3/inputs"
af3_output_path="${exec_directory}/cycle_${i}/af3/outputs"

cd "${work_directory}"
"${biopython_python}" - << EOF
from Bio.PDB import PDBParser, MMCIFIO
structure = PDBParser(QUIET=True).get_structure("protein", "${initial_structure}")
io = MMCIFIO()
io.set_structure(structure)
io.save("${mpnn_input_path}/top_scoring.cif")
EOF

cp "${work_directory}/scripts/biopython_selection.py" "${mpnn_input_path}/"
cp "${work_directory}/scripts/mpnn_to_boltz2_FIXED.sh" "${boltz_input_path}/mpnn_to_boltz2.sh"
cp "${work_directory}/scripts/ipsae.py" "${boltz_output_path}/"
cp "${work_directory}/scripts/surface_area.py" "${boltz_output_path}/"

echo "folders and starting structure ready at $(timestamp)"

echo "step 1: selecting interface residues at $(timestamp)"
cd "${mpnn_input_path}"
"${biopython_python}" biopython_selection.py top_scoring.cif > selection_output.txt
interface_residues=$(grep "^INTERFACE=" selection_output.txt | cut -d'=' -f2)
fixed_residues=$(grep "^FIXED=" selection_output.txt | cut -d'=' -f2)
echo "fixed: ${fixed_residues}"
echo "interface: ${interface_residues}"
cd "${work_directory}"

echo "step 2: running LigandMPNN at $(timestamp)"
sbatch "${work_directory}/submit_mpnn_highmem.sh" \
    "4tool_C3_cycle${i}" \
    "${mpnn_input_path}/top_scoring.cif" \
    "${mpnn_input_path}/../outputs" \
    "${total_sequences}" \
    "${fixed_residues}" \
    "${interface_residues}" \
    "${mpnn_model_type}" \
    "${mpnn_temp}"

sleep 60
desired_lines=$(( 1 + 2 * total_sequences ))
while true; do
    fasta_file=$(find "${mpnn_output_path}" -name "*.fa" 2>/dev/null | head -1)
    if [ -z "${fasta_file:-}" ]; then
        echo "waiting for MPNN output..."
        sleep 30
        continue
    fi
    current_lines=$(wc -l < "${fasta_file}")
    if [ "${current_lines}" -ge "${desired_lines}" ]; then
        echo "MPNN complete at $(timestamp)"
        break
    fi
    sleep 30
done
cp "${fasta_file}" "${trajectory_path}/cycle_${i}/"

echo "step 3: converting FASTA to Boltz-2 YAML (chain-collision-fixed) at $(timestamp)"
bash "${boltz_input_path}/mpnn_to_boltz2.sh" \
    "${fasta_file}" \
    -l "${ligand}:${num_ligands}" \
    -o "${boltz_input_path}/yamls"

echo "step 3b: running Boltz-2 on all sequences at $(timestamp)"
sbatch "${work_directory}/submit_boltz.sh" \
    "${i}_C3" \
    "${boltz_input_path}/yamls" \
    "${boltz_output_path}" \
    "${boltz_samples}"

echo "step 3c: running ESM3 on all sequences at $(timestamp)"
cp "${fasta_file}" "${esm3_input_path}/"
sbatch "${work_directory}/submit_esm3.sh" \
    "${i}_C3" \
    "${esm3_input_path}/$(basename ${fasta_file})" \
    "${esm3_output_path}"

total_boltz_cifs=$(( total_sequences * boltz_samples ))
while true; do
    current_cifs=$(find "${boltz_output_path}" -name "*.cif" 2>/dev/null | wc -l)
    if [ "${current_cifs}" -ge "${total_boltz_cifs}" ]; then
        echo "Boltz-2 complete at $(timestamp)"
        break
    fi
    echo "Boltz-2 running... ${current_cifs}/${total_boltz_cifs}"
    sleep 30
done

while true; do
    if [ -f "${esm3_output_path}/esm3_scores.csv" ]; then
        esm3_lines=$(wc -l < "${esm3_output_path}/esm3_scores.csv")
        if [ "${esm3_lines}" -ge "${total_sequences}" ]; then
            echo "ESM3 complete at $(timestamp)"
            break
        fi
    fi
    echo "ESM3 running..."
    sleep 30
done

echo "step 4: scoring Boltz-2 outputs at $(timestamp)"

"${biopython_python}" - << EOF
import json, csv, glob, os
# FIXED 2026-06-24: previously used Boltz-2's top-level 'iptm', which
# averages across ALL chain pairs (protein + ligand) and is inflated by
# confident haem placement. Now extracts pair_chains_iptm[i][j] for the
# actual protein-protein pair(s), matching Step 7's corrected metric.
# protein_chain_count=3 for this C3 trimer (all three unique pairs averaged).
protein_chain_count = ${protein_chain_count}

def protein_pair_indices(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

pairs = protein_pair_indices(protein_chain_count)

predictions_dir = "${boltz_predictions}"
output_csv = "${boltz_output_path}/scores_file.csv"
rows = []
for seq_dir in glob.glob(f"{predictions_dir}/*/"):
    seq_id = os.path.basename(seq_dir.rstrip('/'))
    json_files = glob.glob(f"{seq_dir}confidence_*.json")
    if not json_files:
        continue
    conf_scores = []
    pp_scores = []
    for jf in json_files:
        with open(jf) as f:
            d = json.load(f)
        conf_scores.append(d.get("confidence_score", 0.0))
        pc = d.get("pair_chains_iptm")
        if pc:
            vals = []
            for i, j in pairs:
                try:
                    v1 = pc[str(i)][str(j)]
                    v2 = pc[str(j)][str(i)]
                    vals.append((v1 + v2) / 2)
                except (KeyError, TypeError):
                    continue
            if vals:
                pp_scores.append(sum(vals) / len(vals))
    mean_conf = sum(conf_scores) / len(conf_scores) if conf_scores else 0.0
    mean_pp_iptm = sum(pp_scores) / len(pp_scores) if pp_scores else 0.0
    rows.append({"id": seq_id, "confidence_score": mean_conf, "iptm": mean_pp_iptm,
                 "calculated_average": (mean_conf + mean_pp_iptm) / 2})
rows.sort(key=lambda x: x["calculated_average"], reverse=True)
with open(output_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["id", "confidence_score", "iptm", "calculated_average"])
    writer.writeheader()
    writer.writerows(rows)
print(f"scored {len(rows)} sequences (using corrected protein-pair iptm)")
EOF

cp "${boltz_output_path}/scores_file.csv" "${trajectory_path}/cycle_${i}/"

echo "step 4b: combining Boltz-2 and ESM3 scores at $(timestamp)"
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

echo "top ${#top_ids[@]} sequences selected for Chai-1 + AF3: ${top_ids[*]}"

echo "step 5: running Chai-1 on top 50 at $(timestamp)"

numeric_ids=()
for id in "${top_ids[@]}"; do
    num=$(echo "${id}" | grep -oP '\d+$')
    [ -n "${num}" ] && numeric_ids+=("${num}")
done

top_ids_pattern=$(printf "| id=%s," "${numeric_ids[@]}")
top_ids_pattern="${top_ids_pattern:2}"
grep -A1 -E "(${top_ids_pattern})" "${fasta_file}" \
    > "${chai_input_path}/top_sequences.fa" 2>/dev/null || true

"${biopython_python}" "${work_directory}/scripts/prep_chai_fasta.py" \
    "${chai_input_path}/top_sequences.fa" \
    "${chai_input_path}/chai_input.fa" \
    "3"

chai_ids_file="${chai_input_path}/chai_ids.txt"
> "${chai_ids_file}"
for id in "${numeric_ids[@]}"; do
    echo "${id}" >> "${chai_ids_file}"
done
n_chai_ids=$(wc -l < "${chai_ids_file}")
echo "submitting Chai-1 array for ${n_chai_ids} ids"

sbatch --array=0-$((n_chai_ids - 1))%5 "${work_directory}/submit_chai_array.sh" \
    "${chai_ids_file}" \
    "${chai_input_path}/chai_input.fa" \
    "${chai_output_path}"

echo "Chai-1 array submitted, waiting for completion at $(timestamp)"
while true; do
    n_chai_jobs=$(squeue -u "$(whoami)" -n chai_array 2>/dev/null | wc -l)
    if [ "${n_chai_jobs}" -le 1 ]; then
        echo "Chai-1 array no longer queued at $(timestamp), verifying completeness..."
        break
    fi
    echo "Chai-1 array running... ${n_chai_jobs} tasks remaining (incl. header)"
    sleep 60
done

missing_chai_ids=()
for id in "${numeric_ids[@]}"; do
    if ! find "${chai_output_path}" -path "*id${id}*combined_scores.csv" 2>/dev/null | grep -q .; then
        missing_chai_ids+=("${id}")
    fi
done
if [ "${#missing_chai_ids[@]}" -gt 0 ]; then
    echo "WARNING: ${#missing_chai_ids[@]} Chai-1 sequences missing output: ${missing_chai_ids[*]}"
    echo "These likely failed (OOM/timeout) and need manual rerun."
    echo "Continuing pipeline with the ${n_chai_ids}-minus-missing sequences that DID complete."
else
    echo "Chai-1 complete at $(timestamp): all ${n_chai_ids} ids present"
fi

echo "step 6: running AF3 (array) on top 50 at $(timestamp)"

"${biopython_python}" "${work_directory}/scripts/fasta_to_af3_nomsa.py" \
    "${chai_input_path}/top_sequences.fa" \
    "${af3_input_path}"

cd "${af3_input_path}"
mkdir -p batches
batch_i=0
batch=0
for f in id*.json; do
    if [ $((batch_i % 5)) -eq 0 ]; then
        batch=$((batch + 1))
        mkdir -p "batches/batch_${batch}"
    fi
    mv "${f}" "batches/batch_${batch}/"
    batch_i=$((batch_i + 1))
done
n_batches=$(ls "${af3_input_path}/batches" | wc -l)
cd "${work_directory}"

echo "submitting AF3 array with ${n_batches} batches"
af3_job_id=$(sbatch --parsable --array=1-${n_batches} "${work_directory}/run_af3_array.sh" \
    "${exec_directory}/cycle_${i}/af3")
echo "AF3 array job ID: ${af3_job_id}"

while true; do
    n_done=$(squeue -j "${af3_job_id}" 2>/dev/null | wc -l)
    if [ "${n_done}" -le 1 ]; then
        echo "AF3 array no longer queued at $(timestamp), verifying completeness..."
        break
    fi
    echo "AF3 array running... ${n_done} tasks remaining (incl. header)"
    sleep 60
done

missing_af3_ids=()
for id in "${numeric_ids[@]}"; do
    if ! find "${af3_output_path}" -path "*id${id}*summary_confidences.json" ! -path "*seed*" 2>/dev/null | grep -q .; then
        missing_af3_ids+=("${id}")
    fi
done
if [ "${#missing_af3_ids[@]}" -gt 0 ]; then
    echo "WARNING: ${#missing_af3_ids[@]} AF3 sequences missing output: ${missing_af3_ids[*]}"
    echo "Check ${work_directory}/logs/ for the corresponding array task logs."
    echo "Continuing pipeline with the sequences that DID complete."
else
    echo "AF3 array complete at $(timestamp): all ${#numeric_ids[@]} ids present"
fi

echo "step 7: compiling CORRECTED scores and selecting next-cycle seeds at $(timestamp)"

"${biopython_python}" - << EOF
import json, csv, glob, os
import numpy as np

protein_chain_count = ${protein_chain_count}

def protein_pair_indices(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

pairs = protein_pair_indices(protein_chain_count)

boltz_predictions = "${boltz_predictions}"
boltz_pp_scores = {}
for seq_dir in glob.glob(f"{boltz_predictions}/*/"):
    seq_id = os.path.basename(seq_dir.rstrip("/")).replace("top_scoring.cif_", "").replace("id", "")
    json_files = glob.glob(f"{seq_dir}confidence_*.json")
    if not json_files:
        continue
    model_pp_vals = []
    for jf in json_files:
        with open(jf) as f:
            d = json.load(f)
        pc = d.get("pair_chains_iptm")
        if not pc:
            continue
        vals = []
        for i, j in pairs:
            try:
                v1 = pc[str(i)][str(j)]
                v2 = pc[str(j)][str(i)]
                vals.append((v1 + v2) / 2)
            except (KeyError, TypeError):
                continue
        if vals:
            model_pp_vals.append(sum(vals) / len(vals))
    if model_pp_vals:
        boltz_pp_scores[seq_id] = max(model_pp_vals)

chai_output_path = "${chai_output_path}"
chai_pp_scores = {}
for npz_path in glob.glob(f"{chai_output_path}/*/scores.model_idx_*.npz"):
    seq_dir = os.path.basename(os.path.dirname(npz_path))
    seq_id = seq_dir.replace("top_scoring.cif_", "").replace("id", "")
    d = np.load(npz_path)
    if "per_chain_pair_iptm" not in d:
        continue
    mat = d["per_chain_pair_iptm"]
    if mat.ndim == 3:
        mat = mat[0]
    vals = []
    for i, j in pairs:
        if i < mat.shape[0] and j < mat.shape[1]:
            vals.append((mat[i, j] + mat[j, i]) / 2)
    if not vals:
        continue
    score = sum(vals) / len(vals)
    if seq_id not in chai_pp_scores or score > chai_pp_scores[seq_id]:
        chai_pp_scores[seq_id] = score

esm3_csv = "${esm3_output_path}/esm3_scores.csv"
esm3_scores = {}
if os.path.exists(esm3_csv):
    with open(esm3_csv) as f:
        for row in csv.DictReader(f):
            seq_id = row["id"].replace("top_scoring.cif_", "").replace("id", "")
            try:
                esm3_scores[seq_id] = float(row["ptm"])
            except (KeyError, ValueError):
                continue

af3_output_path = "${af3_output_path}"
af3_scores = {}
af3_files = glob.glob(f"{af3_output_path}/batch_*/*_summary_confidences.json") + \
            glob.glob(f"{af3_output_path}/batch_*/*/*_summary_confidences.json")
for f in af3_files:
    try:
        with open(f) as fh:
            d = json.load(fh)
        seq_id = os.path.basename(f).split("_summary")[0].replace("id", "")
        vals = []
        for i, j in pairs:
            vals.append(d["chain_pair_iptm"][i][j])
        af3_scores[seq_id] = sum(vals) / len(vals)
    except Exception:
        continue

all_ids = set(chai_pp_scores) | set(af3_scores)
results = []
for seq_id in all_ids:
    c = chai_pp_scores.get(seq_id, 0.0)
    a = af3_scores.get(seq_id, 0.0)
    b = boltz_pp_scores.get(seq_id, 0.0)
    e = esm3_scores.get(seq_id, 0.0)
    results.append({
        "id": seq_id, "boltz_pp": b, "esm3_ptm": e, "chai_pp": c, "af3_pp": a,
        "track1_af3only": a,
        "track2_af3chai": (a + c) / 2,
        "track3_all4": (a + c + b + e) / 4,
    })

ids_sorted_af3 = sorted(results, key=lambda r: -r["af3_pp"])
ids_sorted_chai = sorted(results, key=lambda r: -r["chai_pp"])
af3_rank = {r["id"]: i + 1 for i, r in enumerate(ids_sorted_af3)}
chai_rank = {r["id"]: i + 1 for i, r in enumerate(ids_sorted_chai)}
for r in results:
    r["track4_rank_sum"] = af3_rank[r["id"]] + chai_rank[r["id"]]
    r["track5_chai_only"] = r["chai_pp"]

out_csv = "${trajectory_path}/cycle_${i}/all_scores_CORRECTED.csv"
fieldnames = ["id", "boltz_pp", "esm3_ptm", "chai_pp", "af3_pp",
              "track1_af3only", "track2_af3chai", "track3_all4",
              "track4_rank_sum", "track5_chai_only"]
with open(out_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

print(f"Corrected scores for {len(results)} sequences written to {out_csv}")

seeds_csv = "${trajectory_path}/cycle_${i}/next_cycle_seeds.csv"
seeds_txt = "${trajectory_path}/cycle_${i}/next_cycle_seeds.txt"

track_specs = [
    ("track1_af3only", "Track1_AF3", True),
    ("track2_af3chai", "Track2_AF3Chai_mean", True),
    ("track3_all4", "Track3_All4_mean", True),
    ("track4_rank_sum", "Track4_RankSumConsensus", False),
    ("track5_chai_only", "Track5_Chai_only", True),
]

seed_rows = []
txt_lines = [f"Cycle ${i} (C3) -- next-cycle seed selection per track",
             f"Generated from {len(results)} scored sequences",
             ""]

for track_key, track_name, higher_is_better in track_specs:
    if not results:
        txt_lines.append(f"{track_name}: NO RESULTS FOUND")
        continue
    best = (max if higher_is_better else min)(results, key=lambda r: r[track_key])
    seed_rows.append({
        "track": track_name,
        "winning_id": best["id"],
        "track_score": best[track_key],
        "boltz_pp": best["boltz_pp"],
        "esm3_ptm": best["esm3_ptm"],
        "chai_pp": best["chai_pp"],
        "af3_pp": best["af3_pp"],
    })
    txt_lines.append(
        f"{track_name}: id{best['id']}  "
        f"(track_score={best[track_key]:.4f}, boltz_pp={best['boltz_pp']:.3f}, "
        f"esm3_ptm={best['esm3_ptm']:.3f}, chai_pp={best['chai_pp']:.3f}, "
        f"af3_pp={best['af3_pp']:.3f})"
    )

with open(seeds_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["track", "winning_id", "track_score",
                                            "boltz_pp", "esm3_ptm", "chai_pp", "af3_pp"])
    writer.writeheader()
    writer.writerows(seed_rows)

unique_seed_ids = sorted(set(r["winning_id"] for r in seed_rows), key=lambda x: int(x))
txt_lines.append("")
txt_lines.append(f"Unique seed ids across all tracks ({len(unique_seed_ids)}): "
                  + ", ".join(f"id{x}" for x in unique_seed_ids))

with open(seeds_txt, "w") as f:
    f.write("\n".join(txt_lines) + "\n")

print(f"Per-track seed selection written to {seeds_csv} and {seeds_txt}")
for line in txt_lines:
    print(line)
EOF

echo "C3 cycle ${i} complete at $(timestamp)"
echo "All scores saved to ${trajectory_path}/cycle_${i}/all_scores_CORRECTED.csv"
echo "Next-cycle seeds saved to ${trajectory_path}/cycle_${i}/next_cycle_seeds.csv"
