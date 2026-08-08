#!/usr/bin/env bash
set -u

home_directory="/home/b5ae/mvg2713124.b5ae"
scratch_directory="/scratch/b5ae/mvg2713124.b5ae"
work_directory="${scratch_directory}/cytbx_pipeline"
master_folder="rq2/master_pipeline/RQ2_Q8"
exec_directory="${work_directory}/${master_folder}"
trajectory_path="${exec_directory}/design_trajectory"
biopython_python="${home_directory}/miniconda3/envs/biopython/bin/python"

timestamp() { date +"%F_%T"; }

cofactor="Q8"
top_n_for_chai_af3=50
cofactor_chain_index=3
i=1

boltz_output_path="${exec_directory}/cycle_${i}/boltz/outputs"
boltz_predictions="${boltz_output_path}/boltz_results_input/predictions"
esm3_output_path="${exec_directory}/cycle_${i}/esm3/outputs"
chai_input_path="${exec_directory}/cycle_${i}/chai/inputs"
chai_output_path="${exec_directory}/cycle_${i}/chai/outputs"
af3_input_path="${exec_directory}/cycle_${i}/af3/inputs"
af3_output_path="${exec_directory}/cycle_${i}/af3/outputs"
fasta_file="${exec_directory}/cycle_${i}/boltz/inputs/clashfree_input.fa"

top_ids=()
while IFS=',' read -r id rest; do
    [ "${id}" == "id" ] && continue
    top_ids+=("${id}")
    [ "${#top_ids[@]}" -ge "${top_n_for_chai_af3}" ] && break
done < "${boltz_output_path}/combined_scores.csv"
echo "top ${#top_ids[@]} selected for Chai-1 + AF3: ${top_ids[*]}"

numeric_ids=()
for id in "${top_ids[@]}"; do
    num=$(echo "${id}" | grep -oP '\d+$'); [ -n "${num}" ] && numeric_ids+=("${num}")
done

echo "step 5: Chai-1 on top ${#numeric_ids[@]} at $(timestamp)"
top_ids_pattern=$(printf "| id=%s," "${numeric_ids[@]}"); top_ids_pattern="${top_ids_pattern:2}"
grep -A1 -E "(${top_ids_pattern})" "${fasta_file}" > "${chai_input_path}/top_sequences.fa" 2>/dev/null || true
"${biopython_python}" "${work_directory}/scripts/prep_chai_fasta_q8.py" \
    "${chai_input_path}/top_sequences.fa" \
    "${chai_input_path}/chai_input.fa"
chai_ids_file="${chai_input_path}/chai_ids.txt"; > "${chai_ids_file}"
for id in "${numeric_ids[@]}"; do echo "${id}" >> "${chai_ids_file}"; done
n_chai_ids=$(wc -l < "${chai_ids_file}")
sbatch --array=0-$((n_chai_ids - 1))%5 "${work_directory}/submit_chai_array.sh" \
    "${chai_ids_file}" "${chai_input_path}/chai_input.fa" "${chai_output_path}"

echo "Chai-1 array submitted. Run Step 6 (AF3) separately after Chai-1 completes."
