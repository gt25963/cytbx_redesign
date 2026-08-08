#!/usr/bin/env python3
"""
Convert LigandMPNN fasta output to AF3 JSON inputs (no-MSA).
Usage: python fasta_to_af3_nomsa.py <fasta_path> <output_dir> [n_protein_chains]
  n_protein_chains: oligomeric state (default 2). Ligands = 2 x n_protein (2 HEM per protomer).
"""
import sys, re, json, os, string

fasta_path = sys.argv[1]
output_dir = sys.argv[2]
n_prot = int(sys.argv[3]) if len(sys.argv) > 3 else 2
n_lig = 2 * n_prot                       # 2 HEM per protomer
os.makedirs(output_dir, exist_ok=True)

chain_ids = list(string.ascii_uppercase)  # A, B, C, ...
prot_ids = chain_ids[:n_prot]
lig_ids  = chain_ids[n_prot:n_prot + n_lig]

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

print(f"Found {len(sequences)} sequences; building {n_prot} protein + {n_lig} HEM per design")
for seq_id, seq in sequences:
    seqs = [{"protein": {"id": pid, "sequence": seq, "unpairedMsa": "", "pairedMsa": "", "templates": []}} for pid in prot_ids]
    seqs += [{"ligand": {"id": lid, "ccdCodes": ["HEM"]}} for lid in lig_ids]
    af3_json = {"name": seq_id, "modelSeeds": [1], "sequences": seqs,
                "dialect": "alphafold3", "version": 4}
    with open(os.path.join(output_dir, f"{seq_id}.json"), "w") as f:
        json.dump(af3_json, f)
print(f"Wrote {len(sequences)} AF3 JSON files ({n_prot}-mer) to {output_dir}")
