#!/bin/bash
output_csv="/tmp/burial_results_fmnhemb.csv"
echo "id,cofactor_burial_atoms,total_protein_atoms,reference_burial_atoms,burial_ratio" > "${output_csv}"

for cif in rq2/master_pipeline/RQ2_FMN_hemB/cycle_1/af3/outputs/batch_*/id*/id*_model.cif; do
    id=$(basename "$(dirname "$cif")")
    result=$(python3 /tmp/burial_check.py "$cif" A C FMN B HEM 2>/dev/null)
    cofactor=$(echo "$result" | grep "cofactor_burial_atoms" | awk '{print $2}')
    total=$(echo "$result" | grep "total_protein_atoms" | awk '{print $2}')
    ref=$(echo "$result" | grep "reference_burial_atoms" | awk '{print $2}')
    ratio=$(echo "$result" | grep "burial_ratio_vs_reference" | awk '{print $2}')
    echo "${id},${cofactor},${total},${ref},${ratio}" >> "${output_csv}"
    echo "Processed ${id}: ratio=${ratio}"
done

echo "Done. Results in ${output_csv}"
