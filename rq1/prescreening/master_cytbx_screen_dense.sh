#!/usr/bin/env bash

# CytbX Oligomeric State Screening - DENSE RPXDOCK VARIANT
# Uses pre-existing best dense-sampled RPXDock pose per state
# LigandMPNN -> Boltz-2 + ESM3 -> combined score
# Compares directly against original screen_scores.csv baseline

# This is the 'denser RPXDock sampling' control referenced in Results 1.1.
# Reruns the same screening/scoring pipeline as the original oligomer_screen, but on a single higher-density-sampled RPXDock pose per state (rpxdock_top_n=1 vs the original's top_n=3), to check whether the default sampling density under-sampled the docking space. 
# Produces screen_scores_dense.csv for direct comparison against the original screen_scores.csv.

set -u

# DIRECTORIES
home_directory="/home/b5ae/mvg2713124.b5ae"
scratch_directory="/scratch/b5ae/mvg2713124.b5ae"
work_directory="${scratch_directory}/cytbx_pipeline"
screen_directory="${work_directory}/oligomer_screen_dense"
trajectory_path="${screen_directory}/screen_trajectory"

# PYTHON PATHS
biopython_python="${home_directory}/miniconda3/envs/biopython/bin/python"

# UTILITIES
timestamp() { date +"%F_%T"; }

# INPUT
cytbx_monomer="/home/b5ae/mvg2713124.b5ae/cytbx_redesign/CytbX.pdb"

# PARAMETERS
oligomeric_states=(2 3 4 5)
rpxdock_top_n=1
ligand=HEM
mpnn_model_type="global_label_membrane_mpnn"
mpnn_temp=0.1
total_sequences=20
boltz_samples=2
ipsae_cutoff=15

# FOLDER SETUP
echo "dense screening pipeline starting at $(timestamp)"
mkdir -p "${trajectory_path}"
mkdir -p "${work_directory}/logs"
touch "${trajectory_path}/timestamps.txt"

for state in "${oligomeric_states[@]}"; do
    mkdir -p "${screen_directory}/C${state}/holo"
    for j in $(seq 1 $rpxdock_top_n); do
        mkdir -p "${screen_directory}/C${state}/top${j}/LigandMPNN/inputs"
        mkdir -p "${screen_directory}/C${state}/top${j}/LigandMPNN/outputs"
        mkdir -p "${screen_directory}/C${state}/top${j}/boltz/inputs"
        mkdir -p "${screen_directory}/C${state}/top${j}/boltz/outputs"
        mkdir -p "${screen_directory}/C${state}/top${j}/esm3/inputs"
        mkdir -p "${screen_directory}/C${state}/top${j}/esm3/outputs"
    done
    mkdir -p "${trajectory_path}/C${state}"
done

echo "folders ready at $(timestamp)"

# PART 1: SKIPPED - using pre-existing dense RPXDock poses
# Best pose per state already copied into C${state}/rpxdock/

echo "part 1 skipped - using pre-existing dense RPXDock poses at $(timestamp)"

for state in "${oligomeric_states[@]}"; do
    n_pdbs=$(find "${screen_directory}/C${state}/rpxdock" -name "*.pdb" 2>/dev/null | wc -l)
    echo "C${state}: ${n_pdbs} dense pose(s) found"
done

# PART 2: ADD HAEMS + LIGANDMPNN + BOLTZ-2 + ESM3

for state in "${oligomeric_states[@]}"; do

    echo "processing C${state} at $(timestamp)"
    num_ligands=$(( state * 2 ))

    top_pdbs=$(ls "${screen_directory}/C${state}/rpxdock/"*.pdb 2>/dev/null | \
        grep -v "job0" | sort | head -${rpxdock_top_n})

    j=1
    for pdb in ${top_pdbs}; do

        echo "C${state} top${j}: adding haems and running design"

        holo_pdb="${screen_directory}/C${state}/holo/top${j}_holo.pdb"

        "${biopython_python}" "${work_directory}/scripts/add_haem_to_trimer.py" \
            "${pdb}" \
            "${cytbx_monomer}" \
            "${holo_pdb}"

        mpnn_input="${screen_directory}/C${state}/top${j}/LigandMPNN/inputs"
        mpnn_output="${screen_directory}/C${state}/top${j}/LigandMPNN/outputs/seqs"
        boltz_input="${screen_directory}/C${state}/top${j}/boltz/inputs"
        boltz_output="${screen_directory}/C${state}/top${j}/boltz/outputs"
        esm3_input="${screen_directory}/C${state}/top${j}/esm3/inputs"
        esm3_output="${screen_directory}/C${state}/top${j}/esm3/outputs"

        cp "${holo_pdb}" "${mpnn_input}/top_scoring.pdb"
        cp "${work_directory}/scripts/biopython_selection.py" "${mpnn_input}/"
        cp "${work_directory}/mpnn_to_boltz2.sh" "${boltz_input}/"

        # interface residue selection
        cd "${mpnn_input}"
        "${biopython_python}" biopython_selection.py top_scoring.pdb > selection_output.txt
        interface_residues=$(grep "^INTERFACE=" selection_output.txt | cut -d'=' -f2)
        fixed_residues=$(grep "^FIXED=" selection_output.txt | cut -d'=' -f2)
        cd "${work_directory}"

        # LigandMPNN
        sbatch "${work_directory}/submit_mpnn.sh" \
            "screendense_C${state}_top${j}" \
            "${mpnn_input}/top_scoring.pdb" \
            "${mpnn_input}/../outputs" \
            "${total_sequences}" \
            "${fixed_residues}" \
            "${interface_residues}" \
            "${mpnn_model_type}" \
            "${mpnn_temp}"

        sleep 60

        # wait for FASTA
        desired_lines=$(( 1 + 2 * total_sequences ))
        while true; do
            fasta_file=$(find "${mpnn_output}" -name "*.fa" 2>/dev/null | head -1)
            if [ -z "${fasta_file:-}" ]; then
                echo "waiting for MPNN output C${state} top${j}..."
                sleep 30
                continue
            fi
            current_lines=$(wc -l < "${fasta_file}")
            if [ "${current_lines}" -ge "${desired_lines}" ]; then
                echo "MPNN complete C${state} top${j}"
                break
            fi
            sleep 30
        done

        cp "${fasta_file}" "${trajectory_path}/C${state}/"

        # convert to Boltz-2 YAML
        bash "${boltz_input}/mpnn_to_boltz2.sh" \
            "${fasta_file}" \
            -l "${ligand}:${num_ligands}" \
            -o "${boltz_input}/yamls"

        # Boltz-2
        sbatch "${work_directory}/submit_boltz.sh" \
            "screendense_C${state}_top${j}" \
            "${boltz_input}/yamls" \
            "${boltz_output}" \
            "${boltz_samples}"

        # ESM3
        cp "${fasta_file}" "${esm3_input}/"
        sbatch "${work_directory}/submit_esm3.sh" \
            "screendense_C${state}_top${j}" \
            "${esm3_input}/$(basename ${fasta_file})" \
            "${esm3_output}"

        j=$(( j + 1 ))
    done
done

# wait for all boltz + esm3
echo "waiting for all Boltz-2 and ESM3 runs to complete..."

for state in "${oligomeric_states[@]}"; do
    for j in $(seq 1 $rpxdock_top_n); do
        boltz_output="${screen_directory}/C${state}/top${j}/boltz/outputs"
        esm3_output="${screen_directory}/C${state}/top${j}/esm3/outputs"
        total_cifs=$(( total_sequences * boltz_samples ))

        while true; do
            current_cifs=$(find "${boltz_output}" -name "*.cif" 2>/dev/null | wc -l)
            if [ "${current_cifs}" -ge "${total_cifs}" ]; then
                echo "Boltz-2 C${state} top${j} complete"
                break
            fi
            echo "Boltz-2 C${state} top${j}: ${current_cifs}/${total_cifs}"
            sleep 30
        done

        while true; do
            if [ -f "${esm3_output}/esm3_scores.csv" ]; then
                esm3_lines=$(wc -l < "${esm3_output}/esm3_scores.csv")
                if [ "${esm3_lines}" -ge "${total_sequences}" ]; then
                    echo "ESM3 C${state} top${j} complete"
                    break
                fi
            fi
            echo "waiting for ESM3 C${state} top${j}..."
            sleep 30
        done
    done
done

# PART 3: SCORE AND COMPARE TO ORIGINAL BASELINE

echo "part 3: scoring dense poses at $(timestamp)"

echo "state,top_n,boltz_score,esm3_ptm,combined_score" > "${trajectory_path}/screen_scores_dense.csv"

for state in "${oligomeric_states[@]}"; do
    for j in $(seq 1 $rpxdock_top_n); do
        boltz_output="${screen_directory}/C${state}/top${j}/boltz/outputs"
        esm3_output="${screen_directory}/C${state}/top${j}/esm3/outputs"

        echo "id,confidence_score,protein_iptm,ipSAE,calculated_average" \
            > "${boltz_output}/scores_file.csv"

        cp "${work_directory}/scripts/ipsae.py" "${boltz_output}/"
        cp "${work_directory}/scripts/surface_area.py" "${boltz_output}/"

        for seq_dir in "${boltz_output}/boltz_results_input/predictions/"/*/; do
            seq_id=$(basename "${seq_dir}")
            pae_file=$(find "${seq_dir}" -name "pae_*.npz" 2>/dev/null | head -1)
            cif_file=$(find "${seq_dir}" -name "*.cif" 2>/dev/null | head -1)
            [ -z "${pae_file:-}" ] || [ -z "${cif_file:-}" ] && continue

            "${biopython_python}" "${boltz_output}/ipsae.py" \
                "${pae_file}" "${cif_file}" "${ipsae_cutoff}" "${ipsae_cutoff}"
            "${biopython_python}" "${boltz_output}/surface_area.py" "${cif_file}"

            json_file=$(find "${seq_dir}" -name "confidence_*.json" 2>/dev/null | head -1)
            confidence_score=$(grep -oP '(?<="confidence_score": )[-+]?[0-9]*\.?[0-9]+' "${json_file}")
            protein_iptm=$(grep -oP '(?<="iptm": )[-+]?[0-9]*\.?[0-9]+' "${json_file}")
            ipsae_score=$(awk '$5 == "max" && $6 != "" {sum += $6; count++} \
                END {if (count > 0) printf "%.4f", sum/count}' \
                "${seq_dir}"/*_"${ipsae_cutoff}"_*.txt)
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

        top_combined=$(awk -F',' 'NR==2 {print $4}' "${boltz_output}/combined_scores.csv")
        top_boltz=$(awk -F',' 'NR==2 {print $2}' "${boltz_output}/combined_scores.csv")
        top_esm3=$(awk -F',' 'NR==2 {print $3}' "${boltz_output}/combined_scores.csv")

        echo "C${state},top${j},${top_boltz},${top_esm3},${top_combined}" \
            >> "${trajectory_path}/screen_scores_dense.csv"
    done
done

sort -t, -k5,5 -n -r "${trajectory_path}/screen_scores_dense.csv" \
    -o "${trajectory_path}/screen_scores_dense.csv"

best_state=$(awk -F',' 'NR==2 {print $1}' "${trajectory_path}/screen_scores_dense.csv")
best_score=$(awk -F',' 'NR==2 {print $5}' "${trajectory_path}/screen_scores_dense.csv")

echo "DENSE SCREENING COMPLETE at $(timestamp)"
echo "Best dense oligomeric state: ${best_state} with combined score ${best_score}"
echo "Compare against original baseline in ${work_directory}/screen_scores.csv"

cp "${trajectory_path}/screen_scores_dense.csv" "${work_directory}/"
echo "Dense screen scores saved to ${work_directory}/screen_scores_dense.csv"
