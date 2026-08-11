#!/usr/bin/env bash
set -u

#one-off setup script: builds full folder tree and copies helper scripts into place for fresh C2 test run. COnverts starting scaffold PDB into CIF format for ligmpnn first cycle input
home_directory="/home/b5ae/mvg2713124.b5ae"
scratch_directory="/scratch/b5ae/mvg2713124.b5ae"
work_directory="${scratch_directory}/cytbx_pipeline"
master_folder="CytbX_C2"
exec_directory="${work_directory}/main_pipeline/${master_folder}"
trajectory_path="${exec_directory}/design_trajectory"
biopython_python="${home_directory}/miniconda3/envs/biopython/bin/python"

timestamp() { date +"%F_%T"; }

initial_structure="${work_directory}/prescreening/oligomer_screen/C2/holo/top2_holo.pdb"
number_of_iterations=3

echo "test setup starting at $(timestamp)"
mkdir -p "${exec_directory}"
mkdir -p "${trajectory_path}/initial"
touch "${trajectory_path}/timestamps.txt"
mkdir -p "${work_directory}/logs"

#create full per-tool input/output directory tree for each planned cycle
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

#convert starting monomer PDB into CIF (lignmpnn expected input)
echo "converting starting structure to CIF"
"${biopython_python}" - << PYEOF
from Bio.PDB import PDBParser, MMCIFIO
structure = PDBParser(QUIET=True).get_structure("protein", "${initial_structure}")
io = MMCIFIO()
io.set_structure(structure)
io.save("${exec_directory}/cycle_1/LigandMPNN/inputs/top_scoring.cif")
io.save("${trajectory_path}/initial/top_scoring.cif")
PYEOF

#copy helper scripts into each cycle's working folders so available locally
for i in $(seq 1 $number_of_iterations); do
    cp "${work_directory}/scripts/biopython_selection.py" "${exec_directory}/cycle_${i}/LigandMPNN/inputs/"
    cp "${work_directory}/mpnn_to_boltz2.sh" "${exec_directory}/cycle_${i}/boltz/inputs/"
    cp "${work_directory}/scripts/ipsae.py" "${exec_directory}/cycle_${i}/boltz/outputs/"
    cp "${work_directory}/scripts/surface_area.py" "${exec_directory}/cycle_${i}/boltz/outputs/"
    cp "${work_directory}/scripts/ipsae.py" "${exec_directory}/cycle_${i}/af3/outputs/"
    cp "${work_directory}/scripts/ipsae.py" "${exec_directory}/cycle_${i}/chai/outputs/"
done

echo "setup complete at $(timestamp)"
echo "checking folder structure:"
ls "${exec_directory}"
ls "${exec_directory}/cycle_1/LigandMPNN/inputs/"
