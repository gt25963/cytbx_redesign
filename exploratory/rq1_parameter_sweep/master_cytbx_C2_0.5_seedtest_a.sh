#!/usr/bin/env bash

# CytbX Oligomeric Redesign Pipeline - RQ1 - C2 Dimer
# LigandMPNN -> Boltz-2 (all) + ESM3 (all) -> AF3 + Chai-1 (top N)
# Starting structure: C2 top2 from oligomeric state screening
# This is the 0.5 seedtest_a version

set -u

# DIRECTORIES
home_directory="/home/b5ae/mvg2713124.b5ae"
scratch_directory="/scratch/b5ae/mvg2713124.b5ae"
work_directory="${scratch_directory}/cytbx_pipeline"
master_folder="CytbX_C2_0.5_seedtest_a"
exec_directory="${work_directory}/main_pipeline/${master_folder}"
trajectory_path="${exec_directory}/design_trajectory"
biopython_python="${home_directory}/miniconda3/envs/biopython/bin/python"

# UTILITIES
timestamp() { date +"%F_%T"; }

# STARTING STRUCTURE
initial_structure="/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/prescreening/oligomer_screen/C2/holo/top2_holo.pdb"

# OLIGOMER PARAMETERS
oligomeric_state=2
ligand=HEM
num_ligands=4 # 2 haem per subunit x 2 subunits

# DESIGN PARAMETERS
number_of_iterations=2
total_sequences=50

# LIGANDMPNN
mpnn_model_type="global_label_membrane_mpnn"
mpnn_temp=0.5
mpnn_seed=11111
fixed_residues=""

# STRUCTURE PREDICTION
boltz_samples=2
top_n_for_af3_chai=3
af3_samples=5
chai_samples=5

# SCORING
ipsae_cutoff=15

# FOLDER SETUP
echo "pipeline starting at $(timestamp)"
echo "setting up folder structure for ${number_of_iterations} cycles"

mkdir -p "${exec_directory}"
mkdir -p "${trajectory_path}/initial"
touch "${trajectory_path}/timestamps.txt"
mkdir -p "${work_directory}/logs"

for i in $(seq 1 $number_of_iterations); do
    mkdir -p "${exec_directory}/cycle_${i}/LigandMPNN/inputs"
    mkdir -p "${exec_directory}/cycle_${i}/LigandMPNN/outputs"
    mkdir -p "${exec_directory}/cycle_${i}/boltz/inputs"
    mkdir -p "${exec_directory}/cycle_${i}/boltz/outputs"
    mkdir -p "${exec_directory}/cycle_${i}/esm3/inputs"
    mkdir -p "${exec_directory}/cycle_${i}/esm3/outputs"
    mkdir -p "${exec_directory}/cycle_${i}/af3/inputs"
    mkdir -p "${exec_directory}/cycle_${i}/af3/outputs"
    mkdir -p "${exec_directory}/cycle_${i}/chai/inputs"
    mkdir -p "${exec_directory}/cycle_${i}/chai/outputs"
    mkdir -p "${trajectory_path}/cycle_${i}"
done

# copy starting structure and convert to CIF
cd "${work_directory}"
echo "converting starting structure to CIF"
"${biopython_python}" - << EOF
from Bio.PDB import PDBParser, MMCIFIO
structure = PDBParser(QUIET=True).get_structure("protein", "${initial_structure}")
io = MMCIFIO()
io.set_structure(structure)
io.save("${exec_directory}/cycle_1/LigandMPNN/inputs/top_scoring.cif")
io.save("${trajectory_path}/initial/top_scoring.cif")
EOF

# copy scripts into cycle folders
for i in $(seq 1 $number_of_iterations); do
    cp "${work_directory}/scripts/biopython_selection.py" "${exec_directory}/cycle_${i}/LigandMPNN/inputs/"
    cp "${work_directory}/mpnn_to_boltz2.sh" "${exec_directory}/cycle_${i}/boltz/inputs/"
    cp "${work_directory}/scripts/ipsae.py" "${exec_directory}/cycle_${i}/boltz/outputs/"
    cp "${work_directory}/scripts/surface_area.py" "${exec_directory}/cycle_${i}/boltz/outputs/"
    cp "${work_directory}/scripts/ipsae.py" "${exec_directory}/cycle_${i}/af3/outputs/"
    cp "${work_directory}/scripts/ipsae.py" "${exec_directory}/cycle_${i}/chai/outputs/"
done

echo "folders and starting structure ready at $(timestamp)"

# MAIN ITERATION LOOP
echo "iteration starting at $(timestamp)"
echo "iteration starting at $(timestamp)" >> "${trajectory_path}/timestamps.txt"

for i in $(seq 2 $number_of_iterations); do

    echo "cycle ${i} of ${number_of_iterations} starting at $(timestamp)"
    echo "cycle ${i} starting at $(timestamp)" >> "${trajectory_path}/timestamps.txt"

    mpnn_input_path="${exec_directory}/cycle_${i}/LigandMPNN/inputs"
    mpnn_output_path="${exec_directory}/cycle_${i}/LigandMPNN/outputs/seqs"
    boltz_input_path="${exec_directory}/cycle_${i}/boltz/inputs"
    boltz_output_path="${exec_directory}/cycle_${i}/boltz/outputs"
    boltz_predictions="${boltz_output_path}/boltz_results_input/predictions"
    esm3_input_path="${exec_directory}/cycle_${i}/esm3/inputs"
    esm3_output_path="${exec_directory}/cycle_${i}/esm3/outputs"
    af3_input_path="${exec_directory}/cycle_${i}/af3/inputs"
    af3_output_path="${exec_directory}/cycle_${i}/af3/outputs"
    chai_input_path="${exec_directory}/cycle_${i}/chai/inputs"
    chai_output_path="${exec_directory}/cycle_${i}/chai/outputs"

    # STEP 1: select interface residues 
    echo "step 1: selecting interface residues at $(timestamp)"
    cd "${mpnn_input_path}"
    "${biopython_python}" biopython_selection.py top_scoring.cif > selection_output.txt
    interface_residues=$(grep "^INTERFACE=" selection_output.txt | cut -d'=' -f2)
    fixed_residues_detected=$(grep "^FIXED=" selection_output.txt | cut -d'=' -f2)

    if [ -n "${fixed_residues}" ]; then
        echo "using manually set fixed residues: ${fixed_residues}"
    else
        fixed_residues="${fixed_residues_detected}"
        echo "using auto-detected haem-coordinating residues: ${fixed_residues}"
    fi
    echo "interface residues for design: ${interface_residues}"
    cd "${work_directory}"

    
    # STEP 2: ligMPNN
    echo "step 2: running LigandMPNN cycle ${i} at $(timestamp)"
    sbatch "${work_directory}/submit_mpnn_seeded.sh" \
        "${i}" \
        "${mpnn_input_path}/top_scoring.cif" \
        "${mpnn_input_path}/../outputs" \
        "${total_sequences}" \
        "${fixed_residues}" \
        "${interface_residues}" \
        "${mpnn_model_type}" \
        "${mpnn_temp}" \
        "${mpnn_seed}"

    sleep 60

    desired_lines=$(( 1 + 2 * total_sequences ))
    while true; do
        fasta_file=$(find "${mpnn_output_path}" -name "*.fa" 2>/dev/null | head -1)
        if [ -z "${fasta_file:-}" ]; then
            echo "waiting for MPNN output..."
            sleep 30
            continue
        fi
        current_lines=$(wc -l < "${fasta_file}")
        if [ "${current_lines}" -ge "${desired_lines}" ]; then
            echo "MPNN complete at $(timestamp)"
            break
        fi
        sleep 30
    done

    cp "${fasta_file}" "${trajectory_path}/cycle_${i}/"

    
    # STEP 3: FASTA -> BOLTZ YAML
    echo "step 3: converting FASTA to Boltz-2 YAML at $(timestamp)"
    bash "${boltz_input_path}/mpnn_to_boltz2.sh" \
        "${fasta_file}" \
        -l "${ligand}:${num_ligands}" \
        -o "${boltz_input_path}/yamls"

    
    # STEP 4: BOLTZ-2 + ESM3 parallel
    echo "running Boltz-2 on all sequences at $(timestamp)"
    sbatch "${work_directory}/submit_boltz.sh" \
        "${i}" \
        "${boltz_input_path}/yamls" \
        "${boltz_output_path}" \
        "${boltz_samples}"

    echo "running ESM3 on all sequences at $(timestamp)"
    cp "${fasta_file}" "${esm3_input_path}/"
    sbatch "${work_directory}/submit_esm3.sh" \
        "${i}" \
        "${esm3_input_path}/$(basename ${fasta_file})" \
        "${esm3_output_path}"

    # wait for Boltz-2
    total_boltz_cifs=$(( total_sequences * boltz_samples ))
    while true; do
        current_cifs=$(find "${boltz_output_path}" -name "*.cif" 2>/dev/null | wc -l)
        if [ "${current_cifs}" -ge "${total_boltz_cifs}" ]; then
            echo "Boltz-2 complete at $(timestamp)"
            break
        fi
        echo "Boltz-2 running... ${current_cifs}/${total_boltz_cifs}"
        sleep 30
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

    
    # STEP 5: score Botlz outputs
    echo "step 5: scoring Boltz-2 outputs at $(timestamp)"

    "${biopython_python}" - << EOF
import json, csv, glob, os

predictions_dir = "${boltz_predictions}"
output_csv = "${boltz_output_path}/scores_file.csv"

rows = []
for seq_dir in glob.glob(f"{predictions_dir}/*/"):
    seq_id = os.path.basename(seq_dir.rstrip('/'))
    json_files = glob.glob(f"{seq_dir}confidence_*.json")
    if not json_files:
        continue
    scores = []
    for jf in json_files:
        with open(jf) as f:
            d = json.load(f)
        scores.append(d.get("confidence_score", 0.0))
    mean_conf = sum(scores) / len(scores)
    mean_iptm = sum(json.load(open(jf)).get("iptm", 0.0) for jf in json_files) / len(json_files)
    rows.append({"id": seq_id, "confidence_score": mean_conf, "iptm": mean_iptm,
                 "calculated_average": (mean_conf + mean_iptm) / 2})

rows.sort(key=lambda x: x["calculated_average"], reverse=True)
with open(output_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["id", "confidence_score", "iptm", "calculated_average"])
    writer.writeheader()
    writer.writerows(rows)
print(f"scored {len(rows)} sequences")
EOF

    cp "${boltz_output_path}/scores_file.csv" "${trajectory_path}/cycle_${i}/"

    
    # STEP 6: combine Boltz + ESM3 scores
    echo "step 6: combining scores at $(timestamp)"
    "${biopython_python}" "${work_directory}/scripts/combine_scores.py" \
        "${boltz_output_path}/scores_file.csv" \
        "${esm3_output_path}/esm3_scores.csv" \
        "${boltz_output_path}/combined_scores.csv"

    top_ids=()
    while IFS=',' read -r id rest; do
        [ "${id}" == "id" ] && continue
        top_ids+=("${id}")
        [ "${#top_ids[@]}" -ge "${top_n_for_af3_chai}" ] && break
    done < "${boltz_output_path}/combined_scores.csv"

    echo "top sequences: ${top_ids[*]}"
    cp "${boltz_output_path}/combined_scores.csv" "${trajectory_path}/cycle_${i}/"

    
    # STEP 7: AF3 structure prediction (top N)
    echo "step 7: running AF3 on top ${top_n_for_af3_chai} sequences at $(timestamp)"

    af3_failed=0

    for id in "${top_ids[@]}"; do
        seq=$(grep -A1 "id=${id}" "${fasta_file}" 2>/dev/null | tail -1 | cut -d':' -f1)
        echo -e ">top_scoring,id=${id}\n${seq}" > "${af3_input_path}/${id}.fa"
        "${biopython_python}" "${work_directory}/scripts/mpnn2json_af3.py" \
            "${af3_input_path}/${id}.fa" \
            "${ligand}" \
            "${af3_input_path}"
    done

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
    af3_wait=0

    sleep 60
    if ! ls "${work_directory}/logs/af3_cycle${i}"*.out 2>/dev/null | grep -q "." || \
       grep -r "permission denied" "${work_directory}/logs/af3_cycle${i}"*.out 2>/dev/null | grep -q "weights"; then
        echo "WARNING: AF3 unavailable - skipping AF3 for cycle ${i}"
        af3_failed=1
    fi

    if [ "${af3_failed}" -eq 0 ]; then
        while true; do
            current_af3=$(find "${af3_output_path}" -name "*.cif" 2>/dev/null | wc -l)
            if [ "${current_af3}" -ge "${total_af3_cifs}" ]; then
                echo "AF3 complete at $(timestamp)"
                break
            fi
            af3_wait=$(( af3_wait + 30 ))
            if [ "${af3_wait}" -ge 10800 ]; then
                echo "WARNING: AF3 timed out - skipping AF3 for cycle ${i}"
                af3_failed=1
                break
            fi
            echo "AF3 running... ${current_af3}/${total_af3_cifs}"
            sleep 30
        done
    fi

    if [ "${af3_failed}" -eq 0 ]; then
        cp "${af3_output_path}"/*.cif "${trajectory_path}/cycle_${i}/" 2>/dev/null || true
        echo "AF3 outputs saved to trajectory"
    else
        echo "AF3 skipped - proceeding to Chai-1"
    fi

    
    # STEP 8: Chai predicition
    echo "step 8: running Chai-1 on top ${top_n_for_af3_chai} sequences at $(timestamp)"

    # extract numeric IDs from top_ids for FASTA matching
    numeric_ids=()
    for id in "${top_ids[@]}"; do
        num=$(echo "${id}" | grep -oP '\d+$')
        [ -n "${num}" ] && numeric_ids+=("${num}")
    done

    top_ids_pattern=$(printf "| id=%s," "${numeric_ids[@]}")
    top_ids_pattern="${top_ids_pattern:2}"
    grep -A1 -E "(${top_ids_pattern})" "${fasta_file}" \
        > "${chai_input_path}/top_sequences.fa" 2>/dev/null || true

    "${biopython_python}" "${work_directory}/scripts/prep_chai_fasta.py" \
        "${chai_input_path}/top_sequences.fa" \
        "${chai_input_path}/chai_input.fa" \
        "2"

    sbatch "${work_directory}/submit_chai.sh" \
        "${i}" \
        "${chai_input_path}/chai_input.fa" \
        "${chai_output_path}"

    echo "Chai-1 submitted for cycle ${i} - running in background, pipeline continuing"

    
    # STEP 9: select top structures for next cycle
    echo "step 9: selecting top structure for next cycle at $(timestamp)"

    top_id=$(awk -F',' 'NR==2 {print $1}' "${boltz_output_path}/combined_scores.csv")
    echo "top scoring sequence: ${top_id}"

    top_cif=$(find "${boltz_predictions}" -path "*${top_id}*" -name "*.cif" | head -1)

    if [ "${i}" -lt "${number_of_iterations}" ]; then
        top_seq=$(grep -A1 ", id=${top_id##*_id}," "${fasta_file}" 2>/dev/null | tail -1)
        "${biopython_python}" "${work_directory}/scripts/build_oligomer_cif.py" \
            "${work_directory}/prescreening/oligomer_screen/C2/holo/top2_holo.pdb" \
            "${top_seq}" \
            "${exec_directory}/cycle_$((i+1))/LigandMPNN/inputs/top_scoring.cif"
        echo "top structure built and copied to cycle $((i+1))"
        echo "cycle ${i} complete at $(timestamp)" >> "${trajectory_path}/timestamps.txt"
    else
        top_seq=$(grep -A1 ", id=${top_id##*_id}," "${fasta_file}" 2>/dev/null | tail -1)
        "${biopython_python}" "${work_directory}/scripts/build_oligomer_cif.py" \
            "${work_directory}/prescreening/oligomer_screen/C2/holo/top2_holo.pdb" \
            "${top_seq}" \
            "${trajectory_path}/cycle_${i}/top_scoring_final.cif"
        echo "final cycle complete"
        echo "pipeline complete at $(timestamp)" >> "${trajectory_path}/timestamps.txt"
    fi

    echo "cycle ${i} complete at $(timestamp)"

done

# COMPILE FINAL SCORES
echo "compiling scores across all cycles at $(timestamp)"
echo "cycle,top_id,confidence_score,iptm,combined_score" \
    > "${trajectory_path}/compiled_scores.csv"

for i in $(seq 1 $number_of_iterations); do
    top_line=$(awk -F',' 'NR==2' "${trajectory_path}/cycle_${i}/combined_scores.csv")
    echo "${i},${top_line}" >> "${trajectory_path}/compiled_scores.csv"
done

echo "pipeline complete at $(timestamp)"
echo "final scores at: ${trajectory_path}/compiled_scores.csv"
