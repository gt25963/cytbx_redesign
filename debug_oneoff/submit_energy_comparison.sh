#!/bin/bash
#SBATCH --job-name=energy_comparison
#SBATCH --partition=workq
#SBATCH --account=brics.b5ae
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/energy_comparison_%j.out
conda run -n pyrosetta python /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/diagnostics/energy_comparison/energy_comparison.py
