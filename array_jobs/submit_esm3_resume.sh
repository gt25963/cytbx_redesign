#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=esm3_cycle${1}
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --time=03:00:00 # same as esm3 resumable script but with shorter running time instead of 4 hours
#SBATCH --mem=50GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/esm3_cycle%j.out
#SBATCH --account=brics.b5ae

# Usage: sbatch submit_esm3.sh <cycle_number> <input_fasta> <output_dir>

cycle=$1
input_fasta=$2
output_dir=$3

source /home/b5ae/mvg2713124.b5ae/miniconda3/etc/profile.d/conda.sh
conda activate esm

python /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/scripts/esm3_score_resumable.py \
    "${input_fasta}" \
    "${output_dir}"

echo "ESM3 cycle ${cycle} complete"
