#!/bin/bash
#SBATCH --job-name=AF3
#SBATCH --gpus=1
#SBATCH --time=0:30:00         # Hours:Mins:Secs
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/af3_test/log/%j.out

## Usage: sbatch run_af3.sh <input_path> <output_path>

output_path=$2 
input_path=$1 
mkdir -p $output_path

echo $input_path
weight_path=$SCRATCH/weights
database_path=/projects/b5ae/AF3/database/
af3_container=/projects/b5ae/AF3/af3.sif


singularity exec --nv $af3_container bash -c 'nvidia-smi & nvcc --version'

singularity exec \
    --nv \
    --bind $input_path:/opt/af_input \
    --bind $output_path:/opt/af_output \
    --bind $weight_path:/opt/models \
    --bind $database_path:/opt/public_databases \
    $af3_container \
    python /opt/alphafold3/run_alphafold.py \
        --input_dir=/opt/af_input/ \
        --model_dir=/opt/models \
        --db_dir=/opt/public_databases \
        --output_dir=/opt/af_output