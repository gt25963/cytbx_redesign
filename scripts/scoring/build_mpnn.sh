#!/bin/sh
#SBATCH --job-name=build_mpnn
#SBATCH --output=build_mpnn.out
#SBATCH --gpus=1
#SBATCH --time=12:0:00

# Builds one-off environment - LigMPNN singularity container - once per env setup 
hostname
nvidia-smi

#singularity build boltz.sif boltz_venv.def
singularity build ligandmpnn.sif ligandmpnn.def
