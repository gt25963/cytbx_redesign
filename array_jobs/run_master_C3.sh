#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=master_C3_4tool
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/master_C3_%j.out
#SBATCH --account=brics.b5ae

#EUO: Exit immediately on any error, undefined variable, or failed pipeline step
set -euo pipefail
SCRIPT_DIR="/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/scripts"
#Ensure log directory exists before job starts writing to it
mkdir -p /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs
echo "Launching master_cytbx_4tool_C3.sh at $(date +%F_%T)"
#Run the main C3 4-tool pipeline driver as a background job on the cluster
bash "${SCRIPT_DIR}/master_cytbx_4tool_C3.sh"
echo "master_cytbx_4tool_C3.sh finished at $(date +%F_%T)"
