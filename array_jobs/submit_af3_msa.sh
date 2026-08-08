#!/bin/bash
#SBATCH --partition=workq
#SBATCH --job-name=af3_cycle${1}
#SBATCH --gpus=1
#SBATCH --time=3:00:00
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/%j.out
#SBATCH --account=brics.b5ae

# Usage: sbatch submit_af3.sh <cycle_number> <input_path> <output_path> <json_filename>

cycle=$1
input_path=$2
output_path=$3
json_path=$4

weight_path=/scratch/b5ae/mvg2713124.b5ae/weights
database_path=/projects/b5ae/AF3/database
AF3_path=/projects/b5ae/AF3/alphafold3
af3_container=/projects/b5ae/AF3/AF3_old.sif

singularity exec --nv \
    --bind "${input_path}":/opt/af_input \
    --bind "${output_path}":/opt/af_output \
    --bind "${weight_path}":/opt/models \
    --bind "${database_path}":/opt/public_databases \
    --bind "${AF3_path}":/opt/alphafold3 \
    "${af3_container}" \
    python /opt/alphafold3/run_alphafold.py \
        --json_path=/opt/af_input/"${json_path}" \
        --model_dir=/opt/models \
        --db_dir=/opt/public_databases \
        --output_dir=/opt/af_output
        

echo "AF3 cycle ${cycle} complete for ${json_path}"