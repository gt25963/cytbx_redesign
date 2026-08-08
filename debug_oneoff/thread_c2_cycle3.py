#!/usr/bin/env python3
"""
Thread id11 sequence onto relaxed C2 backbone for cycle 3 input.
"""
from Bio import SeqIO
import re

# Get id11 sequence from cycle 2 combined_scores.csv / fasta
id11_fasta = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/main_pipeline/CytbX_4tool/cycle_2/chai/outputs/top_scoring.cif_id11/combined_scores.csv"

# Read id11 design sequence from the Chai input fasta
chai_input = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/main_pipeline/CytbX_4tool/cycle_2/chai/inputs/chai_input.fa"
id11_seq = None
for record in SeqIO.parse(chai_input, "fasta"):
    if "id11" in record.id or "id=11" in record.id:
        id11_seq = str(record.seq).split("|")[0].strip()
        break

if not id11_seq:
    print("ERROR: Could not find id11 sequence in chai_input.fa")
    exit(1)

print(f"id11 sequence: {id11_seq}")

# Load relaxed backbone
from Bio.PDB import PDBParser, PDBIO

parser = PDBParser()
structure = parser.get_structure("id11_relaxed", "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/main_pipeline/CytbX_4tool/cycle_2/id11_relaxed_membrane.pdb")
model = structure[0]

# Thread sequence onto backbone (keep backbone CA/C/N/O, update sidechains via rotamers)
# For simplicity, use build_oligomer_cif.py logic: substitute sequence
threaded_pdb = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/main_pipeline/CytbX_4tool/cycle_2/id11_relaxed_threaded.pdb"

# Use PyRosetta to repack
from pyrosetta import *
init()

pose = pose_from_file("/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/main_pipeline/CytbX_4tool/cycle_2/id11_relaxed_membrane.pdb")

# Set sequence (for chain A, since C2 is only 2 chains of protein)
# Assuming chain A is the designable chain
for i, aa in enumerate(id11_seq, start=1):
    if i <= pose.total_residue():
        res = pose.residue(i)
        if res.is_polymer():
            mutate_residue(pose, i, aa)

pose.dump_pdb(threaded_pdb)
print(f"Threaded structure saved: {threaded_pdb}")

