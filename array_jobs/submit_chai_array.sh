#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=chai_array
#SBATCH --nodes=1
#SBATCH --time=1:30:00
#SBATCH --mem=256GB
#SBATCH --array=0-44%5
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/chai_array_%A_%a.out
#SBATCH --account=brics.b5ae
# Usage: sbatch submit_chai_array.sh <ids_list_file> <chai_input_fa> <output_base_dir>
ids_list=$1
chai_input_fa=$2
output_base_dir=$3
line_num=$((SLURM_ARRAY_TASK_ID + 1))
target_id=$(sed -n "${line_num}p" "${ids_list}")
if [ -z "${target_id}" ]; then
    echo "ERROR: no id found at line ${line_num} of ${ids_list}"; exit 1
fi
echo "Task ${SLURM_ARRAY_TASK_ID}: processing id=${target_id}"
# find this id's chain1 header line (prefix-agnostic)
start_line=$(grep -nE "name=.*_id${target_id}_chain1$" "${chai_input_fa}" | head -1 | cut -d: -f1)
if [ -z "${start_line}" ]; then
    echo "ERROR: could not find id=${target_id} in ${chai_input_fa}"; exit 1
fi
# find the NEXT design's chain1 header to bound this block (handles any chain count)
next_line=$(awk -v s="${start_line}" 'NR>s && /name=.*_chain1$/ {print NR; exit}' "${chai_input_fa}")
if [ -z "${next_line}" ]; then
    end_line=$(wc -l < "${chai_input_fa}")        # last block: go to EOF
else
    end_line=$((next_line - 1))
fi
echo "Block for id=${target_id}: lines ${start_line}-${end_line} ($((end_line-start_line+1)) lines)"
task_input_dir="$(dirname "${output_base_dir}")/array_inputs"
mkdir -p "${task_input_dir}"
task_fasta="${task_input_dir}/id${target_id}.fa"
sed -n "${start_line},${end_line}p" "${chai_input_fa}" > "${task_fasta}"
echo "protein chains extracted: $(grep -c '^>protein' "${task_fasta}")"
source /home/b5ae/mvg2713124.b5ae/miniconda3/etc/profile.d/conda.sh
conda activate chai
python /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/scripts/batch_aggregate_rank_v2.py \
    "${task_fasta}" "${output_base_dir}"
echo "Task ${SLURM_ARRAY_TASK_ID} (id=${target_id}) complete"
