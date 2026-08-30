#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=screen_dense_master
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=8GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/screen_dense_master_%j.out
#SBATCH --account=brics.b5ae

# run_screen_dense.sh
# same as run_screen.sh but for the denser RPXDock sampling control. 
bash /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/prescreening/master_cytbx_screen_dense.sh
