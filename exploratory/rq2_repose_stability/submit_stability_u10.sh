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

# submit_stability_u10.sh  (RQ2, rq2_repose_stability - exploratory)
# NOTE: Submission wrapper for stability_worker_u10.py, specific to the discarded U10/Q10 quinone naming (superseded by Q8). 
# Not part of the reported pipeline; kept for reference only.
# usage: sbatch --array=0-N%10 submit_stability_u10.sh <ids_file>

ids_file=$1
work=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline

# Each array task picks its own design_id by line number from ids_file, so one submission covers the whole batch
design_id=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${ids_file}")
echo "task ${SLURM_ARRAY_TASK_ID}: U10 id${design_id} at $(date +%T)"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pyrosetta
cd "${work}"

python rq2/scripts/stability_worker_u10.py "${design_id}"
echo "done id${design_id} at $(date +%T)"
