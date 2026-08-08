#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=cytbx_C3_seedtest_c
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=32GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/cytbx_C3_seedtest_c_%j.out
#SBATCH --account=brics.b5ae
bash /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/main_pipeline/master_cytbx_C3_0.3_seedtest_c.sh
