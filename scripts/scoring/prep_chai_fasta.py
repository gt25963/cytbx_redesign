#!/usr/bin/env python3
"""
Convert LigandMPNN FASTA output to Chai-1 input format.
Usage: python prep_chai_fasta.py <input.fasta> <output.fasta> <oligomeric_state>
"""

from pathlib import Path
import re
import sys

def reformat_fasta(input_fasta, output_fasta, oligomeric_state):
    ligand_smiles = "CC1=C(CCC(O)=O)C2=[N]3C1=Cc1c(C)c(C=C)c4C=C5C(C)=C(C=C)C6=[N]5[Fe]3(n14)n1c(=C6)c(C)c(CCC(O)=O)c1=C2"
    
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        lines = [l.strip() for l in infile.readlines()]
        i = 0
        ligand_counter = 1
        while i < len(lines):
            line = lines[i]
            if line == '--' or line == '':
                i += 1
                continue
            if line.startswith(">"):
                match = re.search(r">([^,]+).*id=(\d+)", line)
                if match and i + 1 < len(lines):
                    protein_name = match.group(1).strip()
                    protein_id = match.group(2).strip()
                    seq_line = lines[i + 1]
                    # use : chains if present, otherwise duplicate by oligomeric state
                    if ':' in seq_line:
                        chains = seq_line.split(":")
                    else:
                        chains = [seq_line] * oligomeric_state
                    for chain_idx, chain_seq in enumerate(chains):
                        outfile.write(f">protein|name={protein_name}_id{protein_id}_chain{chain_idx+1}\n")
                        outfile.write(f"{chain_seq}\n")
                        outfile.write(f">ligand|name=heme_b_{ligand_counter}\n{ligand_smiles}\n")
                        ligand_counter += 1
                        outfile.write(f">ligand|name=heme_b_{ligand_counter}\n{ligand_smiles}\n")
                        ligand_counter += 1
                    i += 2
                    continue
                else:
                    print(f"Skipping invalid header: {line}")
            i += 1
    print(f"Reformatted FASTA saved to {output_fasta}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python prep_chai_fasta.py <input.fasta> <output.fasta> <oligomeric_state>")
        sys.exit(1)

    reformat_fasta(sys.argv[1], sys.argv[2], int(sys.argv[3]))