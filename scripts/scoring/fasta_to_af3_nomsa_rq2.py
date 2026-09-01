#!/usr/bin/env python3

# RQ2 AF3 JSON prep: one protein chain (A) + retained HEM (B) + cofactor (C).
# This was used for FMN - there is stale u10 but kept in for the sake of transparency 

import sys, re, json, os

fasta_path, output_dir, cofactor = sys.argv[1], sys.argv[2], sys.argv[3]
if cofactor not in ("FMN", "U10"): ## U10 was before q8 - this script kept because this was used for the FMN
    sys.exit(f"Unknown cofactor {cofactor}; expected FMN or U10")
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
        "name": seq_id,
        "modelSeeds": [1],
        "sequences": [
            {"protein": {"id": "A", "sequence": seq, "unpairedMsa": "", "pairedMsa": "", "templates": []}},
            {"ligand": {"id": "B", "ccdCodes": ["HEM"]}},
            {"ligand": {"id": "C", "ccdCodes": [cofactor]}},
        ],
        "dialect": "alphafold3",
        "version": 4,
    }
    with open(os.path.join(output_dir, f"{seq_id}.json"), "w") as f:
        json.dump(af3, f)
print(f"Wrote {len(sequences)} AF3 JSON files to {output_dir}")
