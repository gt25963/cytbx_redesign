#!/usr/bin/env python3

# ESM3 scoring, resumable. 
# Skips designs whose pdb already exists, and appends each design to the CSV as it completes 

import os, sys, re, csv, torch
from huggingface_hub import login
from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein, GenerationConfig
os.environ['TOKENIZERS_PARALLELISM'] = 'false'
os.environ['HF_DATASETS_CACHE'] = '/scratch/b5ae/mvg2713124.b5ae/.cache/huggingface/'
login(token="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") ## this token has now been deleted so no longer valid - would need to replace with live token before rerunning script

def parse_name(header):
    header = header.lstrip('>')
    base = header.split(',')[0].strip()
    m = re.search(r'id=(\d+)', header)
    return f"{base}_id{m.group(1)}" if m else base

def main(fasta_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    scores_file = os.path.join(output_dir, "esm3_scores.csv")
    # ids already in the CSV 
    done = set()
    if os.path.exists(scores_file):
        with open(scores_file) as f:
            for row in csv.DictReader(f):
                done.add(row["id"])
    write_header = not os.path.exists(scores_file)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Loading ESM3 on {device}")
    model = ESM3.from_pretrained("esm3_sm_open_v1", device=device).float()

    with open(fasta_file) as f:
        arr = [l.strip() for l in f if l.strip()]
    names = [parse_name(arr[i]) for i in range(0, len(arr), 2)]
    seqs  = [arr[i].split(':')[0] for i in range(1, len(arr), 2)]

    csvf = open(scores_file, "a", newline="")
    w = csv.DictWriter(csvf, fieldnames=["id","ptm","plddt"])
    if write_header: w.writeheader(); csvf.flush()

    for seq, name in zip(seqs, names):
        pdb_out = os.path.join(output_dir, f"{name}.pdb")
        if name in done:
            print(f"skip {name} (already done: csv={name in done}, pdb={os.path.exists(pdb_out)})"); continue
        print(f"Processing {name}...")
        protein = ESMProtein(sequence=seq)
        protein = model.generate(protein,
            GenerationConfig(track="structure", num_steps=25, temperature=0.1))
        protein.to_pdb(pdb_out)
        ptm = float(protein.ptm) if protein.ptm is not None else 0.0
        plddt = float(protein.plddt.mean()) if protein.plddt is not None else 0.0
        print(f"{name} — pTM: {ptm:.4f}, mean pLDDT: {plddt:.4f}")
        w.writerow({"id": name, "ptm": f"{ptm:.4f}", "plddt": f"{plddt:.4f}"})
        csvf.flush() ## checkpoint every design
    csvf.close()
    print(f"Done. Scores in {scores_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python esm3_score_resumable.py <input.fasta> <output_dir>"); sys.exit(1)
    if not os.path.isfile(sys.argv[1]):
        print(f"Error: {sys.argv[1]} does not exist."); sys.exit(1)
    main(sys.argv[1], sys.argv[2])
