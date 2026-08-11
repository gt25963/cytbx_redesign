#!/bin/bash
#SBATCH --job-name=fastdesign_C2
#SBATCH --partition=workq
#SBATCH --account=brics.b5ae
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=20:00:00
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/fastdesign_C2_%j.out

#submit script for rosetta approach 3 - applied to c2 dimer, using approach 2 relaxed structure as input. Generate 20 independent decoys
conda run -n pyrosetta python /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/approach3/fastdesign_membrane.py \
    /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/approach2/top2_relaxed_apo.pdb \
    /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/main_pipeline/CytbX_C2_approach3/fastdesign_decoys \
    /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/cytbx_C2_top2_dimer.span \
    20
