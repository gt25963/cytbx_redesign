#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=mpnn_cycle${1}
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --time=2:00:00
#SBATCH --mem=32GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/mpnn_cycle%j.out
#SBATCH --account=brics.b5ae

cycle=$1
input_cif=$2
output_dir=$3
num_sequences=$4
fixed_residues=$5
interface_residues=$6
model_type=$7
temperature=$8
seed=$9

source /home/b5ae/mvg2713124.b5ae/miniconda3/etc/profile.d/conda.sh
conda activate ligandmpnn

cd /scratch/b5ae/mvg2713124.b5ae/LigandMPNN

# same as submit_mpnn.sh, but specificially with seed argument for reproducilbe/ seed-controlled design runs (e.g. seedtest variants in exploratory/)
python run.py \
    --model_type "${model_type}" \
    --pdb_path "${input_cif}" \
    --out_folder "${output_dir}" \
    --fixed_residues "${fixed_residues}" \
    --redesigned_residues "${interface_residues}" \
    --temperature "${temperature}" \
    --seed "${seed}" \
    --pack_side_chains 1 \
    --pack_with_ligand_context 1 \
    --number_of_packs_per_design 1 \
    --global_transmembrane_label 1 \
    --homo_oligomer 1 \
    --batch_size "${num_sequences}"

echo "LigandMPNN cycle ${cycle} complete"
