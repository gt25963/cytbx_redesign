#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=esm3_resume
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --time=4:00:00
#SBATCH --mem=50GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/esm3_resume_%j.out
#SBATCH --account=brics.b5ae
# Usage: sbatch submit_esm3_resumable.sh <cycle_number> <input_fasta> <output_dir>
cycle=$1
input_fasta=$2
output_dir=$3
# activate esm3
source /home/b5ae/mvg2713124.b5ae/miniconda3/etc/profile.d/conda.sh
conda activate esm

#ESM3 scoring - skip any sequences already scored in output directory (not just CSV record, check existing PDB output)
#Use this instead of submit_esm3.sh if previous run timed out partway through the batch
python /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/scripts/esm3_score_resumable.py \
    "${input_fasta}" \
    "${output_dir}"
echo "ESM3 resumable cycle ${cycle} complete"
