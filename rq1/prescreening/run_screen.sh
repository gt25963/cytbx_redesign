#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=cytbx_screen
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=8GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/screen_%j.out
#SBATCH --account=brics.b5ae

source /home/b5ae/mvg2713124.b5ae/miniconda3/etc/profile.d/conda.sh
bash /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/master_cytbx_screen.sh
