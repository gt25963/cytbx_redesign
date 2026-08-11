#!/usr/bin/env bash
#SBATCH --job-name=af3_array
#SBATCH --partition=workq
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --output=af3_array_%a.log

#Exit immediately if an undefined variable is referenced (catch typos early)
set -u

input_fasta="$1"
output_dir="$2"
task_id=${SLURM_ARRAY_TASK_ID}

# Create output directory
mkdir -p "${output_dir}"

# Extract the task_id-th sequence from the fasta
seq_num=$((task_id + 1))
sequence=$(sed -n "$((2 * seq_num))p" "${input_fasta}")
header=$(sed -n "$((2 * seq_num - 1))p" "${input_fasta}")

# Pull design id from header, fall back to task_id if not found
design_id=$(echo "${header}" | grep -oP 'id=\K\d+' || echo "${task_id}")

echo "Running AF3 on ${header} (task ${task_id})"

# Run AF3 on single sequence
cd /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline
python3 -m alphafold.run \
    --fasta_paths="${input_fasta}" \
    --output_dir="${output_dir}/id${design_id}" \
    --model_preset=multimer_v3 \
    --max_sequence_length=2048 \
    --db_preset=full_dbs \
    --job_index="${task_id}"

echo "AF3 complete for id${design_id}"
