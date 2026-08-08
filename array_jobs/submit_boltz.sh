#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=boltz_cycle${1}
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --time=4:00:00
#SBATCH --mem=50GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/boltz_cycle%j.out
#SBATCH --account=brics.b5ae

# Usage: sbatch submit_boltz.sh <cycle_number> <yaml_folder> <output_dir> <diffusion_samples>

cycle=$1
yaml_folder=$2
output_dir=$3
diffusion_samples=$4

# libcuda symlink fix
USER_CUDA_COMPAT="$HOME/boltz_cuda_compat"
REAL_LIBCUDA="/usr/lib64/libcuda.so.565.57.01"
mkdir -p "$USER_CUDA_COMPAT/lib"
if [ ! -f "$USER_CUDA_COMPAT/lib/libcuda.so.1" ]; then
    ln -s "$REAL_LIBCUDA" "$USER_CUDA_COMPAT/lib/libcuda.so.1"
fi

singularity exec --nv \
    --bind "$yaml_folder":/input \
    --bind "$output_dir":/output \
    --bind "$USER_CUDA_COMPAT/lib":/usr/local/cuda/compat \
    /projects/b5ae/containers/boltz_venv.sif \
    boltz predict /input \
        --out_dir /output \
        --accelerator gpu \
        --output_format mmcif \
        --diffusion_samples "${diffusion_samples}" \
        --write_full_pae \
        --no_kernels

echo "Boltz-2 cycle ${cycle} complete"