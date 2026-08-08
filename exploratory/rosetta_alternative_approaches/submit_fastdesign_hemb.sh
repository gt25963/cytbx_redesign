#!/bin/bash
#SBATCH --account=brics.b5ae
#SBATCH --partition=workq
#SBATCH --job-name=fastdes_fmnB
#SBATCH --time=16:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/fastdes_fmnB_%j.out
source /home/b5ae/mvg2713124.b5ae/miniconda3/etc/profile.d/conda.sh
conda activate pyrosetta
cd /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline
python rq2/design/FMN_hemB_pocket/fastdesign_fmn_hemb.py \
  rq2/design/FMN_hemB_pocket/holo_FMN_hemB_prerelax.pdb \
  rq2/design/FMN_hemB_pocket/fastdesign_out \
  rq2/design/FMN_hemB_pocket/cytbx_monomer.span \
  8
echo "fastdesign done"
