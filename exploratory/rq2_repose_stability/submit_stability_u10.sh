#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=stab_u10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=1:00:00
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/stab_u10_%A_%a.out
#SBATCH --account=brics.b5ae

ids_file=$1
work=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline
design_id=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${ids_file}")
echo "task ${SLURM_ARRAY_TASK_ID}: U10 id${design_id} at $(date +%T)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pyrosetta
cd "${work}"
python rq2/scripts/stability_worker_u10.py "${design_id}"
echo "done id${design_id} at $(date +%T)"
