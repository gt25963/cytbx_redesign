#!/bin/bash
# Run ipSAE on all C2 Boltz model_0 structures, collate the 'max' row ipSAE.
set -u
ROOT=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline
PRED=$ROOT/main_pipeline/CytbX_4tool/cycle_1/boltz/outputs/boltz_results_input/predictions
OUT=$ROOT/main_pipeline/CytbX_4tool/cycle_1/ipsae_c2_scores.csv

# Activate conda environment with Biopython/scoring dependencies
source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh
conda activate biopython 2>/dev/null || conda activate base

# Set up output CSV header
echo "id,ipsae_max,ipTM_af,pdockq,pdockq2,LIS" > "$OUT"

# Loop through every candidate's Boltz output directory
for d in "$PRED"/top_scoring.cif_id*/; do
    id=$(basename "$d" | grep -oP 'id\K[0-9]+')
    npz="$d/pae_top_scoring.cif_id${id}_model_0.npz"
    cif="$d/top_scoring.cif_id${id}_model_0.cif"
    # Skip candidates missing required PAE or structure files
    [ -f "$npz" ] && [ -f "$cif" ] || { echo "id$id: missing files" >&2; continue; }
    # Run ipSAE scoring script on this candidate's PAE and structure
    ( cd "$d" && python "$ROOT/scripts/ipsae.py" \
        "pae_top_scoring.cif_id${id}_model_0.npz" \
        "top_scoring.cif_id${id}_model_0.cif" 15 15 >/dev/null 2>&1 )
    txt=$(ls "$d"/*_15_15.txt 2>/dev/null | head -1)
    if [ -n "$txt" ]; then
        # The 'max' row: $5=="max"; cols: ipSAE=$6 ipTM_af=$9 pDockQ=$11 pDockQ2=$12 LIS=$13
        line=$(awk '$5=="max"{print $6","$9","$11","$12","$13; exit}' "$txt")
        echo "id${id},${line}" >> "$OUT"
    else
        echo "id$id: no output" >&2
    fi
done
echo "wrote $OUT"

# Show top 15 candidates ranked by ipSAE score
sort -t, -k2,2 -nr "$OUT" | head -15
