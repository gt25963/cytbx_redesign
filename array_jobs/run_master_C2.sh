#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=master_C2_4tool
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/master_C2_rerun_%j.out
#SBATCH --account=brics.b5ae

#EUO: Exit immediately on any error, undefined variable, or failed pipeline step
set -euo pipefail

#Ensure log directory exists before job starts writing to it
mkdir -p /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs
echo "Launching corrected master_cytbx_4tool.sh (C2 full rerun) at $(date +%F_%T)"

#Run the main C2 4-tool pipeline driver as a background job on the cluster
bash "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/main_pipeline/master_cytbx_4tool.sh"
echo "master_cytbx_4tool.sh finished at $(date +%F_%T)"
