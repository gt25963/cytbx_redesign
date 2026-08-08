#!/bin/bash
#SBATCH --job-name=fastrelax_C3
#SBATCH --partition=workq
#SBATCH --account=brics.b5ae
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/fastrelax_C3_%j.out

conda run -n pyrosetta python /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/approach2/fastrelax_membrane.py \
    /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/prescreening/oligomer_screen/C3/rpxdock/cytbx_C3_CytbX__top3_3.pdb \
    /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/approach2/C3/top3_relaxed_apo.pdb \
    /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/cytbx_C3_top3_trimer.span
