#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=cytbx_C2_0.7
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=32GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/cytbx_C2_0.7_%j.out
#SBATCH --account=brics.b5ae
bash /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/main_pipeline/master_cytbx_C2_0.7.sh
