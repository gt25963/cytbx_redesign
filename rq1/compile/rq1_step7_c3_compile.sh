#!/usr/bin/env bash
set -u

# Compiles Boltz-2/ESM3/Chai-1/AF3 scores for every design in a cycle into one CSV and reports the top design under 5 different track definitions, mirroring the FMN/Q8 compile scripts but using protein-protein interface pairs (0,1)/(0,2)/(1,2) instead of a protein-cofactor pair.

home="/home/b5ae/mvg2713124.b5ae"
work="/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline"
base="${work}/main_pipeline/CytbX_4tool/cycle_2"
traj="${work}/main_pipeline/CytbX_4tool_C3/design_trajectory/cycle_2"
biopython="${home}/miniconda3/envs/biopython/bin/python"
mkdir -p "${traj}"

echo "Step 7: compiling 5 tracks for RQ1 cycle_2, $(date)"

"${biopython}" "${work}/scripts/rq1_step7_c3_compile.py" "${base}" "${traj}"

echo "Step 7 complete, $(date)"
