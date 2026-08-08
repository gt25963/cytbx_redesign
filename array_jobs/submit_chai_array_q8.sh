#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=chai_rq2
#SBATCH --nodes=1
#SBATCH --time=1:30:00
#SBATCH --mem=128GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/chai_rq2_%A_%a.out
#SBATCH --account=brics.b5ae
# Usage: sbatch submit_chai_array_rq2.sh <ids_list> <chai_input_fa> <output_base_dir>
ids_list=$1
chai_input_fa=$2
output_base_dir=$3
line_num=$((SLURM_ARRAY_TASK_ID + 1))
target_id=$(sed -n "${line_num}p" "${ids_list}")
[ -z "${target_id}" ] && { echo "ERROR: no id at line ${line_num} of ${ids_list}"; exit 1; }
echo "Task ${SLURM_ARRAY_TASK_ID}: processing id=${target_id}"
# prefix-agnostic: match any header ending in _id<N>_chain1
start_line=$(grep -nE "name=.*_id${target_id}_chain1$" "${chai_input_fa}" | head -1 | cut -d: -f1)
[ -z "${start_line}" ] && { echo "ERROR: could not find id=${target_id} in ${chai_input_fa}"; exit 1; }
# Q8 record = protein + HEM_B + HEM_C + Q8 = 4 records = 8 lines
end_line=$((start_line + 7))
task_input_dir="${output_base_dir}/array_inputs"
mkdir -p "${task_input_dir}"
task_fasta="${task_input_dir}/id${target_id}.fa"
sed -n "${start_line},${end_line}p" "${chai_input_fa}" > "${task_fasta}"
echo "Extracted lines ${start_line}-${end_line} for id=${target_id}"
source /home/b5ae/mvg2713124.b5ae/miniconda3/etc/profile.d/conda.sh
conda activate chai
python /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/scripts/batch_aggregate_rank_v2.py \
    "${task_fasta}" "${output_base_dir}"
echo "Task ${SLURM_ARRAY_TASK_ID} (id=${target_id}) complete"
