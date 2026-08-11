#!/usr/bin/env bash
#SBATCH --account=brics.b5ae
#SBATCH --partition=workq
#SBATCH --job-name=resume_rq2
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=8GB
#SBATCH --output=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/logs/resume_rq2_%j.out

#catches typos early
set -u
# Usage: sbatch resume_rq2_step4.sh <FMN|U10>
# resumes RQ2 master from Step 4 onward. Asserts Boltz+ESM3 complete (does NOT re-run them).
variant="$1"
home_directory="/home/b5ae/mvg2713124.b5ae"
scratch_directory="/scratch/b5ae/mvg2713124.b5ae"
work_directory="${scratch_directory}/cytbx_pipeline"
biopython_python="${home_directory}/miniconda3/envs/biopython/bin/python"
timestamp() { date +"%F_%T"; }

#Set cofactor-specific paths (FMN or discontinued U10 track)
if [ "${variant}" == "FMN" ]; then
    master_folder="rq2/master_pipeline/RQ2_FMN"; cofactor="FMN"
elif [ "${variant}" == "U10" ]; then
    master_folder="rq2/master_pipeline/RQ2_U10"; cofactor="U10"
else
    echo "Usage: sbatch resume_rq2_step4.sh <FMN|U10>"; exit 1
fi

exec_directory="${work_directory}/${master_folder}"
trajectory_path="${exec_directory}/design_trajectory"
cofactor_chain_index=2
top_n_for_chai_af3=50
i=1
boltz_input_path="${exec_directory}/cycle_${i}/boltz/inputs"
boltz_output_path="${exec_directory}/cycle_${i}/boltz/outputs"
# OVERRIDE: real dir is boltz_results_input, not boltz_results_yamls
boltz_predictions="${boltz_output_path}/boltz_results_input/predictions"
esm3_output_path="${exec_directory}/cycle_${i}/esm3/outputs"
chai_input_path="${exec_directory}/cycle_${i}/chai/inputs"
chai_output_path="${exec_directory}/cycle_${i}/chai/outputs"
af3_input_path="${exec_directory}/cycle_${i}/af3/inputs"
af3_output_path="${exec_directory}/cycle_${i}/af3/outputs"
fasta_file="${esm3_output_path}/../inputs/$(ls ${exec_directory}/cycle_${i}/esm3/inputs/ 2>/dev/null | grep '\.fa$' | head -1)"
mkdir -p "${chai_input_path}" "${chai_output_path}" "${af3_input_path}" "${af3_output_path}" "${trajectory_path}/cycle_${i}"

# ASSERT Boltz complete (do NOT re-run) 
#count confidence JSONs against expected yaml count; abort rather than silently continuing on incomplete data
n_yamls=$(ls "${boltz_input_path}/yamls/"*.yaml 2>/dev/null | wc -l)
n_conf=$(find "${boltz_predictions}" -name "confidence_*.json" 2>/dev/null | wc -l)
echo "[${variant}] Boltz check: ${n_conf} confidence jsons, ${n_yamls} yamls (expect >= ${n_yamls})"
if [ "${n_conf}" -lt "${n_yamls}" ]; then
    echo "ABORT: Boltz incomplete (${n_conf} < ${n_yamls} designs). Not resuming."; exit 1
fi
echo "Boltz complete (${n_conf} jsons for ${n_yamls} designs)."

# ASSERT ESM3 complete (do NOT re-run) 
esm3_csv="${esm3_output_path}/esm3_scores.csv"
if [ ! -f "${esm3_csv}" ]; then echo "ABORT: ${esm3_csv} missing. ESM3 not done."; exit 1; fi
esm3_lines=$(($(wc -l < "${esm3_csv}") - 1))
echo "[${variant}] ESM3 check: ${esm3_lines} designs scored (yamls=${n_yamls})"
if [ "${esm3_lines}" -lt "${n_yamls}" ]; then
    echo "ABORT: ESM3 incomplete (${esm3_lines} < ${n_yamls})."; exit 1
fi
echo "ESM3 complete. fasta_file=${fasta_file}"

echo "[${variant}] entering Step 4 at $(timestamp)"

# Hand off to the master's Step 4 - 7, sourced with our corrected vars 
#export all derived paths/variables so the downstream Step4-7 script can pick up exactly where the resume left off, without re-deriving them itself
export variant home_directory scratch_directory work_directory biopython_python
export exec_directory trajectory_path cofactor cofactor_chain_index top_n_for_chai_af3 i
export boltz_input_path boltz_output_path boltz_predictions esm3_output_path
export chai_input_path chai_output_path af3_input_path af3_output_path fasta_file
bash "${work_directory}/rq2_step4_onward_${variant}.sh"
