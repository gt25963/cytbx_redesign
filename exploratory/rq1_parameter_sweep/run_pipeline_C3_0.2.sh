#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=cytbx_C3_0.2
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=8GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/pipeline_C2_%j.out
#SBATCH --account=brics.b5ae

bash /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/main_pipeline/master_cytbx_C3_0.2.sh
