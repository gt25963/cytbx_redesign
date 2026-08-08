#!/usr/bin/env bash

# CytbX Oligomeric Redesign Pipeline - RQ1 
# RPXDock -> LigandMPNN -> Boltz-2 (all) -> ESM3 (all) -> AF3 -> Chai-1 (top N)

set -euo pipefail

# DIRECTORIES
# home_directory: miniconda/tool installs only
# everything else: scratch (more storage)

home_directory="/home/b5ae/mvg2713124.b5ae"
scratch_directory="/scratch/b5ae/mvg2713124.b5ae"
work_directory="${scratch_directory}/cytbx_pipeline"
master_folder="CytbX_C3"
exec_directory="${work_directory}/${master_folder}"
trajectory_path="${exec_directory}/design_trajectory"

# UTILITIES
timestamp() {
    date +"%F_%T"
}

# STARTING STRUCTURE
# RPXDock output of CytbX oligomer — .cif preferred, .pdb will be auto-converted
initial_structure="rpx_trimer_holo.pdb"

# OLIGOMER PARAMETERS
oligomeric_state=3    # C3 trimer
ligand=HEM
num_ligands=6         # 2 haem per subunit x 3 subunits

# DESIGN PARAMETERS
number_of_iterations=3    # number of full redesign cycles
total_sequences=20        # LigandMPNN sequences generated per cycle

# LIGANDMPNN
mpnn_model_type="global_label_membrane_mpnn"
mpnn_temp=0.1             # increase if sequence diversity is too low
mpnn_bias_for=""          # fill if biasing toward certain residues
mpnn_bias_against=""
fixed_residues=""         # haem-coordinating His residues from geometry script - e.g. A1,A5,..."

# STRUCTURE PREDICTION
boltz_samples=2           # Boltz-2 models per sequence — fast filter
top_n_for_af3_chai=3      # top sequences passed to AF3 + Chai-1
af3_samples=5             # AF3 models per top sequence
chai_samples=5            # Chai-1 models per top sequence

# SCORING
ipsae_cutoff=15           # angstrom cutoff for interface ipSAE


# FOLDER SETUP
echo "pipeline starting at $(timestamp)"
echo "setting up folder structure for ${number_of_iterations} cycles"

# create main run directory
mkdir -p "${exec_directory}"
cd "${exec_directory}"

# create cycle folders for each tool
for i in $(seq 1 $number_of_iterations); do
    mkdir -p cycle_${i}/LigandMPNN/inputs
    mkdir -p cycle_${i}/LigandMPNN/outputs
    mkdir -p cycle_${i}/boltz/inputs
    mkdir -p cycle_${i}/boltz/outputs
    mkdir -p cycle_${i}/esm3/inputs
    mkdir -p cycle_${i}/esm3/outputs
    mkdir -p cycle_${i}/af3/inputs
    mkdir -p cycle_${i}/af3/outputs
    mkdir -p cycle_${i}/chai/inputs
    mkdir -p cycle_${i}/chai/outputs
done

# create trajectory folder for tracking results across cycles
mkdir -p "${trajectory_path}"
touch "${trajectory_path}/timestamps.txt"
for i in $(seq 1 $number_of_iterations); do
    mkdir -p "${trajectory_path}/cycle_${i}"
done
mkdir -p "${trajectory_path}/initial"

# copy starting structure into place
cd "${work_directory}"
if [[ $initial_structure == *.pdb ]]; then
    echo "converting .pdb to .cif"
    source "${home_directory}/miniconda3/etc/profile.d/conda.sh"
    conda activate biopython
    python <<EOF
from Bio.PDB import PDBParser, MMCIFIO
structure = PDBParser(QUIET=True).get_structure("protein", "$initial_structure")
io = MMCIFIO()
io.set_structure(structure)
io.save("initial_structure.cif")
EOF
    cp initial_structure.cif "${trajectory_path}/initial/"
    cp initial_structure.cif "${exec_directory}/cycle_1/LigandMPNN/inputs/top_scoring.cif"
else
    cp "${initial_structure}" "${trajectory_path}/initial/"
    cp "${initial_structure}" "${exec_directory}/cycle_1/LigandMPNN/inputs/top_scoring.cif"
fi

echo "folders and starting structure ready at $(timestamp)"

# copy pipeline scripts into each cycle folder
for i in $(seq 1 $number_of_iterations); do
    cp "${work_directory}/scripts/biopython_selection.py" "${exec_directory}/cycle_${i}/LigandMPNN/inputs/"
    cp "${work_directory}/mpnn_to_boltz2.sh" "${exec_directory}/cycle_${i}/boltz/inputs/"
    cp "${work_directory}/scripts/ipsae.py" "${exec_directory}/cycle_${i}/boltz/outputs/"
    cp "${work_directory}/scripts/surface_area.py" "${exec_directory}/cycle_${i}/boltz/outputs/"
    cp "${work_directory}/scripts/ipsae.py" "${exec_directory}/cycle_${i}/af3/outputs/"
    cp "${work_directory}/scripts/ipsae.py" "${exec_directory}/cycle_${i}/chai/outputs/"
done

# create logs folder
mkdir -p "${work_directory}/logs"

echo "scripts copied into cycle folders at $(timestamp)"

# =======================
# MAIN ITERATION LOOP
# =======================

echo "iteration starting at $(timestamp)"
echo "iteration starting at $(timestamp)" >> "${trajectory_path}/timestamps.txt"

for i in $(seq 1 $number_of_iterations); do

    echo "cycle ${i} of ${number_of_iterations} starting at $(timestamp)"
    echo "cycle ${i} starting at $(timestamp)" >> "${trajectory_path}/timestamps.txt"

    mpnn_input_path="${exec_directory}/cycle_${i}/LigandMPNN/inputs"
    mpnn_output_path="${exec_directory}/cycle_${i}/LigandMPNN/outputs/seqs"
    boltz_input_path="${exec_directory}/cycle_${i}/boltz/inputs"
    boltz_output_path="${exec_directory}/cycle_${i}/boltz/outputs"
    af3_input_path="${exec_directory}/cycle_${i}/af3/inputs"
    af3_output_path="${exec_directory}/cycle_${i}/af3/outputs"
    chai_input_path="${exec_directory}/cycle_${i}/chai/inputs"
    chai_output_path="${exec_directory}/cycle_${i}/chai/outputs"
    esm3_input_path="${exec_directory}/cycle_${i}/esm3/inputs"
    esm3_output_path="${exec_directory}/cycle_${i}/esm3/outputs"

    # STEP 1: Select Interface Residues
    echo "step 1: selecting interface residues at $(timestamp)"
    cd "${mpnn_input_path}"
    source "${home_directory}/miniconda3/etc/profile.d/conda.sh"
    conda activate biopython
    python biopython_selection.py top_scoring.cif > selection_output.txt
    interface_residues=$(grep "^INTERFACE=" selection_output.txt | cut -d'=' -f2)
    fixed_residues_detected=$(grep "^FIXED=" selection_output.txt | cut -d'=' -f2)

    # use manually set fixed_residues if provided, otherwise use detected haem His
    if [ -n "${fixed_residues}" ]; then
        echo "using manually set fixed residues: ${fixed_residues}"
    else
        fixed_residues="${fixed_residues_detected}"
        echo "using auto-detected haem-coordinating residues: ${fixed_residues}"
    fi
    echo "interface residues for design: ${interface_residues}"

    # STEP 2: LigandMPNN
    echo "step 2: running LigandMPNN cycle ${i} at $(timestamp)"
    sbatch "${work_directory}/submit_mpnn.sh" \
        "${i}" \
        "${mpnn_input_path}/top_scoring.cif" \
        "${mpnn_input_path}/../outputs" \
        "${total_sequences}" \
        "${fixed_residues}" \
        "${interface_residues}" \
        "${mpnn_model_type}" \
        "${mpnn_temp}"

    # wait for MPNN output directory
    while [ ! -d "${mpnn_output_path}" ]; do
        echo "waiting for MPNN output directory..."
        sleep 30
    done

    # wait for complete FASTA
    desired_lines=$(( 1 + 2 * total_sequences ))
    while true; do
        fasta_file=$(find "${mpnn_output_path}" -name "*.fa" | head -1)
        if [ -z "${fasta_file:-}" ]; then
            echo "waiting for FASTA file..."
            sleep 15
            continue
        fi
        current_lines=$(wc -l < "${fasta_file}")
        if [ "${current_lines}" -ge "${desired_lines}" ]; then
            echo "MPNN complete — FASTA ready at $(timestamp)"
            break
        else
            echo "FASTA not complete yet (${current_lines}/${desired_lines} lines)..."
            sleep 30
        fi
    done

    cp "${fasta_file}" "${trajectory_path}/cycle_${i}/"
    echo "FASTA saved to trajectory"

    # STEP 3: Convert Fasta to Boltz-2 YAML
    echo "step 3: converting FASTA to Boltz-2 YAML inputs at $(timestamp)"
    bash "${boltz_input_path}/mpnn_to_boltz2.sh" \
        "${fasta_file}" \
        -l "${ligand}:${num_ligands}" \
        -o "${boltz_input_path}/yamls"

    # STEP 4: Boltz-2 Structure Prediction (all sequences) 
    echo "step 4: running Boltz-2 on all ${total_sequences} sequences at $(timestamp)"
    sbatch "${work_directory}/submit_boltz.sh" \
        "${i}" \
        "${boltz_input_path}/yamls" \
        "${boltz_output_path}" \
        "${boltz_samples}"

    # STEP 4b: ESM3 structure prediction (all seqs) — parallel with Boltz-2
    echo "step 4b: running ESM3 on all ${total_sequences} sequences at $(timestamp)"
    cp "${fasta_file}" "${esm3_input_path}/"
    sbatch "${work_directory}/submit_esm3.sh" \
        "${i}" \
        "${esm3_input_path}/$(basename ${fasta_file})" \
        "${esm3_output_path}"

    # wait for Boltz-2
    while [ ! -d "${boltz_output_path}" ]; do
        echo "waiting for Boltz-2 output directory..."
        sleep 15
    done

    total_boltz_cifs=$(( total_sequences * boltz_samples ))
    while true; do
        current_cifs=$(find "${boltz_output_path}" -name "*.cif" | wc -l)
        if [ "${current_cifs}" -ge "${total_boltz_cifs}" ]; then
            echo "Boltz-2 complete — ${current_cifs} cifs at $(timestamp)"
            break
        else
            echo "Boltz-2 running... ${current_cifs}/${total_boltz_cifs} cifs"
            sleep 30
        fi
    done

    # wait for ESM3
    while true; do
        if [ -f "${esm3_output_path}/esm3_scores.csv" ]; then
            esm3_lines=$(wc -l < "${esm3_output_path}/esm3_scores.csv")
            if [ "${esm3_lines}" -ge "${total_sequences}" ]; then
                echo "ESM3 complete at $(timestamp)"
                break
            fi
        fi
        echo "ESM3 running..."
        sleep 30
    done

    # STEP 5: Score Boltz-2 Outputs
    echo "step 5: scoring Boltz-2 outputs at $(timestamp)"
    echo "id,confidence_score,protein_iptm,ipSAE,calculated_average" \
        > "${boltz_output_path}/scores_file.csv"

    conda activate biopython

    for seq_dir in "${boltz_output_path}"/*/; do
        seq_id=$(basename "${seq_dir}")

        pae_file=$(find "${seq_dir}" -name "pae_*.npz" | head -1)
        cif_file=$(find "${seq_dir}" -name "*.cif" | head -1)

        if [ -z "${pae_file:-}" ] || [ -z "${cif_file:-}" ]; then
            echo "skipping ${seq_id} — missing pae or cif"
            continue
        fi

        python "${boltz_output_path}/ipsae.py" \
            "${pae_file}" "${cif_file}" "${ipsae_cutoff}" "${ipsae_cutoff}"

        python "${boltz_output_path}/surface_area.py" "${cif_file}"

        json_file=$(find "${seq_dir}" -name "confidence_*.json" | head -1)
        confidence_score=$(grep -oP '(?<="confidence_score": )[-+]?[0-9]*\.?[0-9]+' "${json_file}")
        protein_iptm=$(grep -oP '(?<="protein_iptm": )[-+]?[0-9]*\.?[0-9]+' "${json_file}")
        ipsae_score=$(awk '$5 == "max" && $6 != "" {sum += $6; count++} \
            END {if (count > 0) printf "%.4f", sum/count}' \
            "${seq_dir}"/*_"${ipsae_cutoff}"_*.txt)
        calculated_average=$(printf "%.4f" \
            "$(echo "scale=6; ($confidence_score + $protein_iptm + (1.5 * $ipsae_score)) / 3" | bc -l)")

        echo "${seq_id},${confidence_score},${protein_iptm},${ipsae_score},${calculated_average}" \
            >> "${boltz_output_path}/scores_file.csv"
    done

    # sort by score descending
    sort -t, -k5,5 -n -r "${boltz_output_path}/scores_file.csv" \
        -o "${boltz_output_path}/scores_file.csv"
    cp "${boltz_output_path}/scores_file.csv" "${trajectory_path}/cycle_${i}/"
    echo "Boltz-2 scoring complete at $(timestamp)"

    # STEP 6: COmbine Boltz-2 and ESM3 Scores - select top N 
    echo "step 6: combining Boltz-2 and ESM3 scores at $(timestamp)"

    source "${home_directory}/miniconda3/etc/profile.d/conda.sh"
    conda activate biopython

    python3 "${work_directory}/scripts/combine_scores.py" \
            "${boltz_output_path}/scores_file.csv" \
            "${esm3_output_path}/esm3_scores.csv" \
            "${boltz_output_path}/combined_scores.csv"

    # select top N from combined scores
    top_ids=()
    while IFS=',' read -r id rest; do
        [ "${id}" == "id" ] && continue
        top_ids+=("${id}")
        [ "${#top_ids[@]}" -ge "${top_n_for_af3_chai}" ] && break
    done < "${boltz_output_path}/combined_scores.csv"

    echo "top sequences selected: ${top_ids[*]}"
    cp "${boltz_output_path}/combined_scores.csv" "${trajectory_path}/cycle_${i}/"

    # STEP 7: AF3 Structure Prediction (Top N)
    echo "step 7: preparing and running AF3 on top ${top_n_for_af3_chai} sequences at $(timestamp)"

    for id in "${top_ids[@]}"; do
        # extract sequence for this id from fasta
        seq=$(grep -A1 "id=${id}" "${fasta_file}" | tail -1 | cut -d':' -f1)
        echo -e ">top_scoring,id=${id}\n${seq}" > "${af3_input_path}/${id}.fa"
        python "${work_directory}/scripts/mpnn2json_af3.py" \
            "${af3_input_path}/${id}.fa" \
            "${ligand}" \
            "${af3_input_path}"
    done

    # submit one AF3 job per json
    mkdir -p "${af3_output_path}/log"
    for json_file in "${af3_input_path}"/*.json; do
        json_name=$(basename "${json_file}")
        sbatch "${work_directory}/submit_af3.sh" \
            "${i}" \
            "${af3_input_path}" \
            "${af3_output_path}" \
            "${json_name}"
    done

    total_af3_cifs=$(( top_n_for_af3_chai * af3_samples ))
    while true; do
        current_af3=$(find "${af3_output_path}" -name "*.cif" | wc -l)
        if [ "${current_af3}" -ge "${total_af3_cifs}" ]; then
            echo "AF3 complete — ${current_af3} cifs at $(timestamp)"
            break
        else
            echo "AF3 running... ${current_af3}/${total_af3_cifs} cifs"
            sleep 30
        fi
    done

    cp "${af3_output_path}"/*.cif "${trajectory_path}/cycle_${i}/" 2>/dev/null || true

    # STEP 8: Chai-1 Prediction (Top N)
    echo "step 8: running Chai-1 on top ${top_n_for_af3_chai} sequences at $(timestamp)"

    # prepare chai input fasta — top N sequences only
    top_ids_pattern=$(printf "|id=%s" "${top_ids[@]}")
    top_ids_pattern="${top_ids_pattern:1}"
    grep -A1 -E "(${top_ids_pattern})" "${fasta_file}" \
        > "${chai_input_path}/top_sequences.fa" || true

    python "${work_directory}/scripts/prep_chai_fasta.py" \
        "${chai_input_path}/top_sequences.fa" \
        "${chai_input_path}/chai_input.fa"

    sbatch "${work_directory}/submit_chai.sh" \
        "${i}" \
        "${chai_input_path}/chai_input.fa" \
        "${chai_output_path}"

    total_chai_cifs=$(( top_n_for_af3_chai * chai_samples ))
    while true; do
        current_chai=$(find "${chai_output_path}" -name "*.cif" | wc -l)
        if [ "${current_chai}" -ge "${total_chai_cifs}" ]; then
            echo "Chai-1 complete — ${current_chai} cifs at $(timestamp)"
            break
        else
            echo "Chai-1 running... ${current_chai}/${total_chai_cifs} cifs"
            sleep 30
        fi
    done

    cp "${chai_output_path}"/*.cif "${trajectory_path}/cycle_${i}/" 2>/dev/null || true

    # STEP 9: Select Top Scoring Structure for Next Cycle 
    echo "step 9: selecting top structure for next cycle at $(timestamp)"

    # top scoring id from Boltz-2 scores (already sorted)
    top_id=$(awk -F',' 'NR==2 {print $1}' "${boltz_output_path}/combined_scores.csv")
    echo "top scoring sequence for cycle ${i}: ${top_id}"

    # find its cif and copy as input for next cycle
    top_cif=$(find "${boltz_output_path}" -path "*${top_id}*" -name "*.cif" | head -1)

    if [ $i -lt $number_of_iterations ]; then
        cp "${top_cif}" "${exec_directory}/cycle_$((i+1))/LigandMPNN/inputs/top_scoring.cif"
        echo "top structure copied to cycle $((i+1)) inputs"
        echo "cycle ${i} complete at $(timestamp)" >> "${trajectory_path}/timestamps.txt"
    else
        cp "${top_cif}" "${trajectory_path}/cycle_${i}/top_scoring_final.cif"
        echo "final cycle complete — top structure saved to trajectory"
        echo "pipeline complete at $(timestamp)" >> "${trajectory_path}/timestamps.txt"
    fi

    echo "cycle ${i} complete at $(timestamp)"

done

# COMPILE FINAL SCORES
echo "compiling scores across all cycles at $(timestamp)"
echo "cycle,top_id,confidence_score,protein_iptm,ipSAE,mean" \
    > "${trajectory_path}/compiled_scores.csv"

for i in $(seq 1 $number_of_iterations); do
    top_line=$(awk -F',' 'NR==2' "${trajectory_path}/cycle_${i}/combined_scores.csv")
    echo "${i},${top_line}" >> "${trajectory_path}/compiled_scores.csv"
done

echo "pipeline complete at $(timestamp)"
echo "final scores at: ${trajectory_path}/compiled_scores.csv"