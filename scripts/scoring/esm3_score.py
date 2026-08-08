#!/usr/bin/env python3
"""
ESM3 structure prediction and scoring for CytbX designed sequences.
Usage: python esm3_score.py <input.fasta> <output_dir>
"""

import os
import sys
import re
import torch
import csv
from huggingface_hub import login
from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein, GenerationConfig

os.environ['TOKENIZERS_PARALLELISM'] = 'false'
os.environ['HF_DATASETS_CACHE'] = '/scratch/b5ae/mvg2713124.b5ae/.cache/huggingface/'

login(token="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")

def parse_name(header):
    header = header.lstrip('>')
    base = header.split(',')[0].strip()
    id_match = re.search(r'id=(\d+)', header)
    if id_match:
        return f"{base}_id{id_match.group(1)}"
    return base

def main(fasta_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    scores_file = os.path.join(output_dir, "esm3_scores.csv")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Loading ESM3 on {device}")
    model = ESM3.from_pretrained("esm3_sm_open_v1", device=device)
    model = model.float()

    with open(fasta_file, 'r') as f:
        seqs_array = [line.strip() for line in f if line.strip()]

    names = [parse_name(seqs_array[i]) for i in range(0, len(seqs_array), 2)]
    seqs = [seqs_array[i].split(':')[0] for i in range(1, len(seqs_array), 2)]

    results = []
    for seq, name in zip(seqs, names):
        print(f"Processing {name}...")
        pdb_out = os.path.join(output_dir, f"{name}.pdb")

        protein = ESMProtein(sequence=seq)
        protein = model.generate(
            protein,
            GenerationConfig(track="structure", num_steps=25, temperature=0.1)
        )
        protein.to_pdb(pdb_out)

        ptm = float(protein.ptm) if protein.ptm is not None else 0.0
        plddt = float(protein.plddt.mean()) if protein.plddt is not None else 0.0

        print(f"{name} — pTM: {ptm:.4f}, mean pLDDT: {plddt:.4f}")
        results.append({"id": name, "ptm": ptm, "plddt": plddt})

    results.sort(key=lambda x: x["ptm"], reverse=True)
    with open(scores_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["id", "ptm", "plddt"])
        writer.writeheader()
        writer.writerows(results)

    print(f"Scores saved to {scores_file}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python esm3_score.py <input.fasta> <output_dir>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_dir = sys.argv[2]

    if not os.path.isfile(fasta_file):
        print(f"Error: {fasta_file} does not exist.")
        sys.exit(1)

    main(fasta_file, output_dir)
