#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=relax_fmnB
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mem=16GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/relax_fmnB_%j.out
#SBATCH --account=brics.b5ae
source /home/b5ae/mvg2713124.b5ae/miniconda3/etc/profile.d/conda.sh
conda activate pyrosetta
cd /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline
python rq2/design/FMN_hemB_pocket/relax_fmn_hemB.py
echo "relax done"
