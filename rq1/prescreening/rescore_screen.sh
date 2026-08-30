#!/usr/bin/env bash
# Rescore existing screening outputs and generate screen_scores_updated.csv

# Rescores the original C2-C5 oligomer screening outputs (Table 1), reading the "iptm" field instead of "protein_iptm" - the latter is always 0.0 for these multi-chain ligand-containing structures.

set -u

home_directory="/home/b5ae/mvg2713124.b5ae"
work_directory="/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline"
screen_directory="${work_directory}/prescreening/oligomer_screen"
trajectory_path="${screen_directory}/screen_trajectory"
biopython_python="${home_directory}/miniconda3/envs/biopython/bin/python"

oligomeric_states=(2 3 4 5)
rpxdock_top_n=3
ipsae_cutoff=15

echo "state,top_n,boltz_score,esm3_ptm,combined_score" > "${trajectory_path}/screen_scores_updated.csv"

for state in "${oligomeric_states[@]}"; do
    for j in $(seq 1 $rpxdock_top_n); do
        boltz_output="${screen_directory}/C${state}/top${j}/boltz/outputs"
        esm3_output="${screen_directory}/C${state}/top${j}/esm3/outputs"
        predictions_dir="${boltz_output}/boltz_results_input/predictions"

        echo "scoring C${state} top${j}..."

        echo "id,confidence_score,protein_iptm,ipSAE,calculated_average" \
            > "${boltz_output}/scores_file.csv"

        cp "${work_directory}/scripts/ipsae.py" "${boltz_output}/"
        cp "${work_directory}/scripts/surface_area.py" "${boltz_output}/"

        for seq_dir in "${predictions_dir}"/*/; do
            seq_id=$(basename "${seq_dir}")
            pae_file=$(find "${seq_dir}" -name "pae_*.npz" 2>/dev/null | head -1)
            cif_file=$(find "${seq_dir}" -name "*.cif" 2>/dev/null | head -1)
            [ -z "${pae_file:-}" ] || [ -z "${cif_file:-}" ] && continue

            "${biopython_python}" "${boltz_output}/ipsae.py" \
                "${pae_file}" "${cif_file}" "${ipsae_cutoff}" "${ipsae_cutoff}"
            "${biopython_python}" "${boltz_output}/surface_area.py" "${cif_file}"

            json_file=$(find "${seq_dir}" -name "confidence_*.json" 2>/dev/null | head -1)
            confidence_score=$(grep -oP '(?<="confidence_score": )[-+]?[0-9]*\.?[0-9]+' "${json_file}")
            ## corrected field: "iptm", not "protein_iptm" (always 0.0 for these structures)
            protein_iptm=$(grep -oP '(?<="iptm": )[-+]?[0-9]*\.?[0-9]+' "${json_file}")
            ipsae_score=$(awk '$5 == "max" && $6 != "" {sum += $6; count++} \
                END {if (count > 0) printf "%.4f", sum/count}' \
                "${seq_dir}"/*_"${ipsae_cutoff}"_*.txt 2>/dev/null)
            ## fall back to 0.0 if ipSAE produced no rows (e.g. no residues within cutoff)
            ipsae_score=${ipsae_score:-0.0}
            calculated_average=$(printf "%.4f" \
                "$(echo "scale=6; ($confidence_score + $protein_iptm + (1.5 * $ipsae_score)) / 3" | bc -l)")
            echo "${seq_id},${confidence_score},${protein_iptm},${ipsae_score},${calculated_average}" \
                >> "${boltz_output}/scores_file.csv"
        done

        sort -t, -k5,5 -n -r "${boltz_output}/scores_file.csv" \
            -o "${boltz_output}/scores_file.csv"

        "${biopython_python}" "${work_directory}/scripts/combine_scores.py" \
            "${boltz_output}/scores_file.csv" \
            "${esm3_output}/esm3_scores.csv" \
            "${boltz_output}/combined_scores.csv"
       
        ## take just the top-ranked row's scores as this pose's representative score
        top_combined=$(awk -F',' 'NR==2 {print $4}' "${boltz_output}/combined_scores.csv")
        top_boltz=$(awk -F',' 'NR==2 {print $2}' "${boltz_output}/combined_scores.csv")
        top_esm3=$(awk -F',' 'NR==2 {print $3}' "${boltz_output}/combined_scores.csv")

        echo "C${state},top${j},${top_boltz},${top_esm3},${top_combined}" \
            >> "${trajectory_path}/screen_scores_updated.csv"
    done
done

sort -t, -k5,5 -n -r "${trajectory_path}/screen_scores_updated.csv" \
    -o "${trajectory_path}/screen_scores_updated.csv"

best_state=$(awk -F',' 'NR==2 {print $1}' "${trajectory_path}/screen_scores_updated.csv")
best_top=$(awk -F',' 'NR==2 {print $2}' "${trajectory_path}/screen_scores_updated.csv")
best_score=$(awk -F',' 'NR==2 {print $5}' "${trajectory_path}/screen_scores_updated.csv")

echo "SCREENING COMPLETE"
echo "Best oligomeric state: ${best_state} (${best_top}) with combined score ${best_score}"

cp "${trajectory_path}/screen_scores_updated.csv" "${work_directory}/"
echo "Screen scores saved to ${work_directory}/screen_scores_updated.csv"
