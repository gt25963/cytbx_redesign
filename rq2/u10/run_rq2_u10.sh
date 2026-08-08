#!/bin/bash
#SBATCH --account=brics.b5ae
#SBATCH --partition=workq
#SBATCH --job-name=master_RQ2_U10
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=8GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/master_RQ2_U10_%j.out
bash /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/master_pipeline/RQ2_U10/master_cytbx_rq2_u10.sh
