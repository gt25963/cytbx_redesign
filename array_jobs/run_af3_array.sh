#!/bin/bash
#SBATCH --job-name=AF3_array
#SBATCH --gpus=1
#SBATCH --time=0:30:00
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/af3_array_%A_%a.out
## Usage: sbatch --array=1-20 run_af3_array.sh <lineage_base_dir>
lineage_base=$1
input_path="${lineage_base}/inputs/batches/batch_${SLURM_ARRAY_TASK_ID}"
output_path="${lineage_base}/outputs/batch_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${output_path}"
weight_path=$SCRATCH/weights
database_path=/projects/b5ae/AF3/database/
af3_container=/projects/b5ae/AF3/af3.sif
singularity exec \
    --nv \
    --bind ${input_path}:/opt/af_input \
    --bind ${output_path}:/opt/af_output \
    --bind ${weight_path}:/opt/models \
    --bind ${database_path}:/opt/public_databases \
    ${af3_container} \
    python /opt/alphafold3/run_alphafold.py \
        --input_dir=/opt/af_input/ \
        --model_dir=/opt/models \
        --db_dir=/opt/public_databases \
        --output_dir=/opt/af_output
