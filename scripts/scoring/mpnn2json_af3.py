#!/usr/bin/env python3

import json 
import sys
import re

def fasta_to_json(fasta_path):
    sequences = []
    current_id = None
    current_seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Save previous sequence if exists
                if current_id is not None:
                    sequences.append({
                        "id": current_id,
                        "sequence": "".join(current_seq)
                    })

                #current_id = line[1:].strip().split(',')[1]
                #current_id = re.search(r'id=(\d+)', line).group(1)
                current_id = line[1:].strip().split(',')[1].split('=')[1].strip()
                current_seq = []
            else:
                current_seq.append(line)

        # Add the last sequence
        if current_id is not None:
            sequences.append({
                "id": current_id,
                "sequence": "".join(current_seq)
            })

    return sequences

def create_af3_fasta(sequences, ligands_list, output_path=""):
    for seq in sequences:
        af3 = {
            "name" : "seq"+seq['id'],
            "modelSeeds": [1],
            "sequences" : [
                    {"protein" : { 
                        "id" : "A",
                        "sequence" : seq['sequence']}}
                    ,
                {   
                    "ligand" : {
                        "id" : "B",
                        "ccdCodes" : ligands_list
                    }}
            ],
            "dialect" : "alphafold3",
            "version": 4
        }

        with open(output_path+'/'+'seq'+seq['id']+'.json', 'w+') as f:
            json.dump(af3, f)

def main():
    fasta_path=sys.argv[1]
    ligands_input=sys.argv[2]
    outpath=sys.argv[3]

    # Parse comma-separated ligands into a list
    ligands_list = [ligand.strip() for ligand in ligands_input.split(',')]

    sequences = fasta_to_json(fasta_path)
    print(f"Found {len(sequences)} sequences")
    print(f"Ligands: {ligands_list}")
    create_af3_fasta(sequences, ligands_list, outpath)

if __name__ == '__main__':
    main()
