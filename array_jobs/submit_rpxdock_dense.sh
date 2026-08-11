#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=rpxdock_C${3}
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=32GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/rpxdock_%j.out
#SBATCH --account=brics.b5ae

input_monomer=$1
output_dir=$2
olig_state=$3

source /home/b5ae/mvg2713124.b5ae/miniconda3/etc/profile.d/conda.sh
conda activate rpxdock

export CXX=/usr/bin/g++-12
export CC=/usr/bin/gcc-12

hscore_dir=/home/b5ae/mvg2713124.b5ae/miniconda3/envs/rpxdock/lib/python3.12/site-packages/rpxdock/data/hscore

mkdir -p "${output_dir}"
cd "${output_dir}"

# denser sampling control (finer cart/ori resolution, keeping top 100 poses instead of 10) tested against the default-density baseline; did not improve pose quality and was not adopted (Table 1, methods in report)
python -mrpxdock \
    --architecture "C${olig_state}" \
    --inputs1 "${input_monomer}" \
    --output_prefix "cytbx_C${olig_state}" \
    --dump_pdbs \
    --hscore_data_dir "${hscore_dir}" \
    --hscore_files small_ilv_h \
    --recenter_input \
    --cart_resl 5 \
    --ori_resl 15 \
    --grid_resolution_cart_angstroms 0.5 \
    --grid_resolution_ori_degrees 0.5 \
    --nout_top 100

echo "RPXDock C${olig_state} docking complete"
