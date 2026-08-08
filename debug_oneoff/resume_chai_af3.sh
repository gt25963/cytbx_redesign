#!/usr/bin/env bash
#SBATCH --account=brics.b5ae
#SBATCH --partition=workq
#SBATCH --job-name=resume_rq1
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=8GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/resume_%j.out
set -u
# Usage: sbatch resume_chai_af3.sh <C2|C3>
variant="$1"

home_directory="/home/b5ae/mvg2713124.b5ae"
scratch_directory="/scratch/b5ae/mvg2713124.b5ae"
work_directory="${scratch_directory}/cytbx_pipeline"
biopython_python="${home_directory}/miniconda3/envs/biopython/bin/python"
timestamp() { date +"%F_%T"; }

if [ "${variant}" == "C2" ]; then
    master_folder="CytbX_4tool"; protein_chain_count=2; chai_oligo=2
elif [ "${variant}" == "C3" ]; then
    master_folder="CytbX_4tool_C3"; protein_chain_count=3; chai_oligo=3
else
    echo "Usage: sbatch resume_chai_af3.sh <C2|C3>"; exit 1
fi

i=1
exec_directory="${work_directory}/main_pipeline/${master_folder}"
trajectory_path="${exec_directory}/design_trajectory"
boltz_output_path="${exec_directory}/cycle_${i}/boltz/outputs"
boltz_predictions="${boltz_output_path}/boltz_results_input/predictions"
esm3_output_path="${exec_directory}/cycle_${i}/esm3/outputs"
chai_input_path="${exec_directory}/cycle_${i}/chai/inputs"
chai_output_path="${exec_directory}/cycle_${i}/chai/outputs"
af3_input_path="${exec_directory}/cycle_${i}/af3/inputs"
af3_output_path="${exec_directory}/cycle_${i}/af3/outputs"

# full 50 numeric ids (for verification + AF3), and the missing subset (for Chai)
mapfile -t numeric_ids < "${chai_input_path}/chai_ids.txt"
missing_file="${chai_input_path}/chai_ids_missing.txt"
n_missing=$(wc -l < "${missing_file}")
echo "[${variant}] resuming at $(timestamp): ${n_missing} missing Chai ids of ${#numeric_ids[@]} total"

# ===== Step 5 (resume): Chai-1 on MISSING ids only =====
if [ "${n_missing}" -gt 0 ]; then
    sbatch --array=0-$((n_missing - 1))%5 "${work_directory}/submit_chai_array.sh" \
        "${missing_file}" \
        "${chai_input_path}/chai_input.fa" \
        "${chai_output_path}"
    echo "Chai-1 resume array submitted at $(timestamp)"
    while true; do
        n_chai_jobs=$(squeue -u "$(whoami)" -n chai_array 2>/dev/null | wc -l)
        if [ "${n_chai_jobs}" -le 1 ]; then
            echo "Chai-1 array no longer queued at $(timestamp), verifying..."
            break
        fi
        echo "Chai-1 array running... ${n_chai_jobs} remaining (incl. header)"
        sleep 60
    done
fi

missing_chai_ids=()
for id in "${numeric_ids[@]}"; do
    if ! find "${chai_output_path}" -path "*id${id}*combined_scores.csv" 2>/dev/null | grep -q .; then
        missing_chai_ids+=("${id}")
    fi
done
if [ "${#missing_chai_ids[@]}" -gt 0 ]; then
    echo "WARNING: ${#missing_chai_ids[@]} Chai-1 still missing: ${missing_chai_ids[*]}"
else
    echo "Chai-1 complete at $(timestamp): all ${#numeric_ids[@]} ids present"
fi

# ===== Step 6: AF3 on ALL 50 (neither master reached AF3) =====
echo "step 6: AF3 array on top ${#numeric_ids[@]} at $(timestamp)"
"${biopython_python}" "${work_directory}/scripts/fasta_to_af3_nomsa.py" \
    "${chai_input_path}/top_sequences.fa" \
    "${af3_input_path}"
cd "${af3_input_path}"
mkdir -p batches
batch_i=0; batch=0
for f in id*.json; do
    if [ $((batch_i % 5)) -eq 0 ]; then
        batch=$((batch + 1)); mkdir -p "batches/batch_${batch}"
    fi
    mv "${f}" "batches/batch_${batch}/"
    batch_i=$((batch_i + 1))
done
n_batches=$(ls "${af3_input_path}/batches" | wc -l)
cd "${work_directory}"

af3_job_id=$(sbatch --parsable --array=1-${n_batches} "${work_directory}/run_af3_array.sh" \
    "${exec_directory}/cycle_${i}/af3")
echo "AF3 array job ID: ${af3_job_id}"
while true; do
    n_done=$(squeue -j "${af3_job_id}" 2>/dev/null | wc -l)
    if [ "${n_done}" -le 1 ]; then
        echo "AF3 array no longer queued at $(timestamp), verifying..."
        break
    fi
    echo "AF3 array running... ${n_done} remaining (incl. header)"
    sleep 60
done

missing_af3_ids=()
for id in "${numeric_ids[@]}"; do
    if ! find "${af3_output_path}" -path "*id${id}*summary_confidences.json" ! -path "*seed*" 2>/dev/null | grep -q .; then
        missing_af3_ids+=("${id}")
    fi
done
if [ "${#missing_af3_ids[@]}" -gt 0 ]; then
    echo "WARNING: ${#missing_af3_ids[@]} AF3 missing: ${missing_af3_ids[*]}"
else
    echo "AF3 complete at $(timestamp): all ${#numeric_ids[@]} ids present"
fi

echo "[${variant}] Chai + AF3 complete at $(timestamp)."
echo "Run Step 7 compile separately once both states are done:"
echo "  bash step7_compile.sh ${variant}"
