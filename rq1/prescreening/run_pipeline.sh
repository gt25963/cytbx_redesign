#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=cytbx_pipeline
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=8GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/pipeline_%j.out
#SBATCH --account=brics.b5ae

# run_pipeline.sh
# Wraps master_cytbx_rq1.sh in a sbatch job. 
# Hence frees up terminal space while running RQ1 on compute node.
source /home/b5ae/mvg2713124.b5ae/miniconda3/etc/profile.d/conda.sh
bash /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/master_cytbx_rq1.sh

echo "Pipeline submitted for completion"
