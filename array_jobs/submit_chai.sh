#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=chai_cycle${1}
#SBATCH --nodes=1
#SBATCH --time=3:00:00
#SBATCH --mem=200GB
#SBATCH --gres=gpu:1
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/chai_cycle%j.out
#SBATCH --account=brics.b5ae

# Usage: sbatch submit_chai.sh <cycle_number> <input_fasta> <output_dir>

cycle=$1
input_fasta=$2
output_dir=$3

source /home/b5ae/mvg2713124.b5ae/miniconda3/etc/profile.d/conda.sh
conda activate chai

python /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/scripts/batch_aggregate_rank_v2.py \
    "${input_fasta}" \
    "${output_dir}"

echo "Chai-1 cycle ${cycle} complete"