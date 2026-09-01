#!/usr/bin/env python3

#Q8 AF3 JSON prep: protein (A) + BOTH hemes (B,C by CCD HEM) + Q8 (D by CCD UQ8). - basically same as the fasta_to_af3_nomsa_rq2.py but for the q8 only
#Q8 = ubiquinone-8, PDB CCD code UQ8 (C49H74O4, 56 PDB entries).
#Chain order: A protein, B Hem1, C Hem2, D UQ8 -> [0][3] = protein<->Q8.

import sys, re, json, os
fasta_path, output_dir = sys.argv[1], sys.argv[2]
os.makedirs(output_dir, exist_ok=True)
sequences = []
with open(fasta_path) as f:
    header = None
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            header = line
        elif header is not None:
            m = re.search(r"id=(\d+)", header)
            if m:
                sequences.append((f"id{m.group(1)}", line.split(":")[0]))
            header = None
print(f"Found {len(sequences)} sequences")
for seq_id, seq in sequences:
    af3 = {
        "name": seq_id, "modelSeeds": [1],
        "sequences": [
            {"protein": {"id": "A", "sequence": seq, "unpairedMsa": "", "pairedMsa": "", "templates": []}},
            {"ligand": {"id": "B", "ccdCodes": ["HEM"]}},
            {"ligand": {"id": "C", "ccdCodes": ["HEM"]}},
            {"ligand": {"id": "D", "ccdCodes": ["UQ8"]}},
        ],
        "dialect": "alphafold3", "version": 4,
    }
    with open(os.path.join(output_dir, f"{seq_id}.json"), "w") as f:
        json.dump(af3, f)
print(f"Wrote {len(sequences)} AF3 JSON files to {output_dir}")
