#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=cytbx_C2_approach2
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=32GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/cytbx_C2_approach2_%j.out
#SBATCH --account=brics.b5ae
bash /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/approach2/master_cytbx_C2_approach2.sh
