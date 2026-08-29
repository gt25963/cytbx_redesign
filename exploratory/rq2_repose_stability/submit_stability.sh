#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=stability
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=1:00:00
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/stability_%A_%a.out
#SBATCH --account=brics.b5ae

# submit_stability.sh  (RQ2, rq2_repose_stability - exploratory)
# NOTE: Submission wrapper for stability_worker.py, part of the earlier, discontinued cofactor stability-scoring approach. 
# Superseded by the LigandMPNN-based redesign pipeline used for the final reported results.
# usage: sbatch --array=0-N%10 submit_stability.sh <COFACTOR> <ids_file>

cofactor=$1
ids_file=$2
work=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline

# Each array task picks its own design_id by line number from ids_file, so one submission covers the whole batch (%10 caps 10 concurrent tasks)
design_id=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${ids_file}")
echo "task ${SLURM_ARRAY_TASK_ID}: ${cofactor} id${design_id} at $(date +%T)"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pyrosetta
cd "${work}"

python rq2/scripts/stability_worker.py "${cofactor}" "${design_id}"
echo "done id${design_id} at $(date +%T)"
