#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=mpnn_rq2
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --time=4:00:00
#SBATCH --mem=128GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/mpnn_rq2_%j.out
#SBATCH --account=brics.b5ae

cycle=$1
input_pdb=$2
output_dir=$3
num_sequences=$4
fixed_residues=$5
redesigned_residues=$6
model_type=$7
temperature=$8

source /home/b5ae/mvg2713124.b5ae/miniconda3/etc/profile.d/conda.sh
conda activate ligandmpnn
cd /scratch/b5ae/mvg2713124.b5ae/LigandMPNN

# RQ2 (FMN/Q8) redesign: same as RQ1 ligandmpnn, but holds fixed_residues (both native haem-coordinating histidines) constant while redesigning the cofactor pocket, and omits --homo_oligomer since RQ2 designs are monomers
python run.py \
    --model_type "${model_type}" \
    --pdb_path "${input_pdb}" \
    --out_folder "${output_dir}" \
    --fixed_residues "${fixed_residues}" \
    --redesigned_residues "${redesigned_residues}" \
    --temperature "${temperature}" \
    --pack_side_chains 1 \
    --pack_with_ligand_context 1 \
    --number_of_packs_per_design 1 \
    --global_transmembrane_label 1 \
    --batch_size "${num_sequences}"

echo "LigandMPNN RQ2 cycle ${cycle} complete"
