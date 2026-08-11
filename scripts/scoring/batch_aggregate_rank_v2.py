from pathlib import Path
import numpy as np
import torch
from chai_lab.chai1 import run_inference
from Bio import SeqIO
import csv
import sys
import re

# paths passed in as arguments from master script
INPUT_FASTA = Path(sys.argv[1])
OUTPUT_BASE_DIR = Path(sys.argv[2])
AGGREGATE_RANK_FILE = OUTPUT_BASE_DIR / "highest_aggregate_scores.csv"

def parse_fasta(fasta_file):
    """
    Parse FASTA into groups of chains belonging to the same sequence.
    Groups are identified by the base ID (e.g. top_scoring.cif_id2) stripping the _chain suffix.
    Returns a list of dicts: {name, chains: [{name, sequence}], ligands: [{name, smiles}]}
    """
    # Walk through every record, grouping chains/ligands that belong to the same design
    groups = {}
    group_order = []
    current_chain = None

    for record in SeqIO.parse(fasta_file, "fasta"):
        description = record.description
        sequence = str(record.seq)

        if "protein|name=" in description:
            chain_name = description.split("name=")[1]
            # strip _chainN suffix to get group key
            base_name = re.sub(r'_chain\d+$', '', chain_name)
            if base_name not in groups:
                groups[base_name] = {"name": base_name, "chains": [], "ligands": []}
                group_order.append(base_name)
            groups[base_name]["chains"].append({"name": chain_name, "sequence": sequence})
            current_chain = base_name

        elif "ligand|name=" in description and current_chain:
            #ligand records are attached to whichever protein group came before them
            ligand_name = description.split("name=")[1]
            groups[current_chain]["ligands"].append({"name": ligand_name, "smiles": sequence})

    return [groups[k] for k in group_order]


def clean_and_expand_data(data):
    #conver numpy score arrays into plain python values/lists for CSV 
    if isinstance(data, np.ndarray):
        if data.ndim == 0:
            return data.item()
        elif data.ndim == 1:
            return list(data)
        elif data.ndim == 2:
            return data.flatten().tolist()
    return data


def extract_scores_to_dict(npz_file):
    #pull key chai confidence metrics out of a single prediction's npz output
    if not npz_file.exists():
        print(f"File {npz_file} not found.")
        return {}

    data = np.load(npz_file)
    scores = {"file": npz_file.stem}

    metrics_to_extract = ["aggregate_score", "iptm", "ptm"]
    for metric in metrics_to_extract:
        if metric in data:
            scores[metric] = clean_and_expand_data(data[metric])

    return scores


protein_data = parse_fasta(INPUT_FASTA)

#run chai structure prediction separately for each design ( group of chains + ligands)
for group in protein_data:
    # rebuild a standalone fasta containing just this design's chains and ligands
    fasta_content = ""
    for chain in group["chains"]:
        fasta_content += f">protein|name={chain['name']}\n{chain['sequence']}\n"
    for ligand in group["ligands"]:
        fasta_content += f">ligand|name={ligand['name']}\n{ligand['smiles']}\n"

    fasta_path = Path(f"./inputs/{group['name']}.fasta")
    output_dir = OUTPUT_BASE_DIR / group['name']
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_path.write_text(fasta_content)

    #run chai interface directly on this single design
    candidates = run_inference(
        fasta_file=fasta_path,
        output_dir=output_dir,
        num_trunk_recycles=3,
        num_diffn_timesteps=200,
        seed=42,
        device=torch.device("cuda" if torch.cuda.is_available() else "cpu"),
        use_esm_embeddings=True,
    )

    #collect scores from every model npz chai produced for this design
    npz_files = list(output_dir.glob("*.npz"))
    combined_scores_file = output_dir / "combined_scores.csv"

    with combined_scores_file.open(mode="w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=["file", "aggregate_score", "iptm", "ptm"])
        writer.writeheader()
        for npz_file in npz_files:
            scores = extract_scores_to_dict(npz_file)
            if scores:
                writer.writerow(scores)

    print(f"Processed {group['name']} - scores saved to {combined_scores_file}")
    #clean the temporary per-desing fasta once predicition is done
    try:
        fasta_path.unlink()
    except FileNotFoundError:
        pass

# pick out the single best scoring model for each design, across all designs, and rank them into one summary file
all_combined_scores = []

for group in protein_data:
    output_dir = OUTPUT_BASE_DIR / group["name"]
    combined_scores_file = output_dir / "combined_scores.csv"

    if combined_scores_file.exists():
        with combined_scores_file.open(mode="r") as file:
            reader = csv.DictReader(file)
            highest_aggregate_row = None
            for row in reader:
                aggregate_score = float(row["aggregate_score"].strip("[]"))
                if highest_aggregate_row is None or aggregate_score > float(highest_aggregate_row["aggregate_score"].strip("[]")):
                    highest_aggregate_row = row

            if highest_aggregate_row:
                highest_aggregate_row["directory"] = group["name"]
                all_combined_scores.append(highest_aggregate_row)

#write final ranked summary, sorted by aggregate score, highest first. 
# The aggegarte score is chai's overall summary score which is inflated by cofactor predicition confidence. 
# The protein-protein/protein-cofactor score used for candidate selection is extracted separately through per_chain_pair_iptm, not this file
with AGGREGATE_RANK_FILE.open(mode="w", newline="") as file:
    writer = csv.DictWriter(file, fieldnames=["directory", "file", "aggregate_score", "iptm", "ptm"])
    writer.writeheader()
    sorted_scores = sorted(all_combined_scores, key=lambda x: float(x["aggregate_score"].strip("[]")), reverse=True)
    writer.writerows(sorted_scores)

print(f"Highest aggregate scores saved to {AGGREGATE_RANK_FILE}")
