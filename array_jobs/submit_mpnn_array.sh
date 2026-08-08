#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=mpnn_array
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --time=4:00:00
#SBATCH --mem=128GB
#SBATCH --array=0-3
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/mpnn_array_%A_%a.out
#SBATCH --account=brics.b5ae

# Usage: sbatch submit_mpnn_array.sh <input_cif> <output_base_dir> <seqs_per_task> <fixed_residues> <interface_residues> <model_type> <temperature>
input_cif=$1
output_base_dir=$2
seqs_per_task=$3
fixed_residues=$4
interface_residues=$5
model_type=$6
temperature=$7

task_output_dir="${output_base_dir}/batch_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${task_output_dir}"

task_seed=$((SLURM_ARRAY_TASK_ID * 100000 + RANDOM))

source /home/b5ae/mvg2713124.b5ae/miniconda3/etc/profile.d/conda.sh
conda activate ligandmpnn
cd /scratch/b5ae/mvg2713124.b5ae/LigandMPNN

python run.py \
    --model_type "${model_type}" \
    --pdb_path "${input_cif}" \
    --out_folder "${task_output_dir}" \
    --fixed_residues "${fixed_residues}" \
    --redesigned_residues "${interface_residues}" \
    --temperature "${temperature}" \
    --seed "${task_seed}" \
    --pack_side_chains 1 \
    --pack_with_ligand_context 1 \
    --number_of_packs_per_design 1 \
    --global_transmembrane_label 1 \
    --homo_oligomer 1 \
    --batch_size "${seqs_per_task}"

echo "Task ${SLURM_ARRAY_TASK_ID} complete, seed ${task_seed}"
