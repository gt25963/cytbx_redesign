#!/usr/bin/env python3
"""
FastRelax id11 with membrane scoring to improve C2 axis alignment.
First convert CIF to PDB and remove ligands.
"""

from pyrosetta import *
from pyrosetta.rosetta.protocols.relax import FastRelax
import numpy as np
import os

# Convert CIF to PDB first (remove ligands)
from Bio.PDB import PDBParser, PDBIO, Select

class ProteinSelect(Select):
    def accept_residue(self, residue):
        # Only accept standard polymer residues (chains A, B)
        if residue.get_parent().get_id() in ['A', 'B']:
            return 1
        return 0

cif_file = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/main_pipeline/CytbX_4tool/cycle_2/chai/outputs/top_scoring.cif_id11/pred.model_idx_0.cif"
pdb_file = "/tmp/id11_protein_only.pdb"

# Load CIF with Bio.PDB
from Bio.PDB import MMCIFParser
parser = MMCIFParser()
structure = parser.get_structure("id11", cif_file)

# Write PDB with only protein chains
io = PDBIO()
io.set_structure(structure)
io.save(pdb_file, ProteinSelect())

print(f"Converted CIF to PDB (protein only): {pdb_file}")

# Now load with PyRosetta
init(extra_options="-mp:lipids:has_pore false")
pose = pose_from_file(pdb_file)

print(">>> FastRelax id11 with membrane scoring")
print(f"Input pose: {pose.total_residue()} residues")

# Get initial axis Z-component
ca_a = []
ca_b = []
for i in range(1, pose.total_residue() + 1):
    res = pose.residue(i)
    if not res.is_polymer():
        continue
    if res.chain() == 1:
        ca_a.append([res.atom("CA").xyz().x, res.atom("CA").xyz().y, res.atom("CA").xyz().z])
    elif res.chain() == 2:
        ca_b.append([res.atom("CA").xyz().x, res.atom("CA").xyz().y, res.atom("CA").xyz().z])

ca_a = np.array(ca_a)
ca_b = np.array(ca_b)
centroid_a = ca_a.mean(axis=0)
centroid_b = ca_b.mean(axis=0)
axis = centroid_b - centroid_a
axis_norm = axis / np.linalg.norm(axis)
z_initial = abs(axis_norm[2])

print(f"Initial axis Z-component: {z_initial:.4f}")

# Setup scoring with membrane context
scorefxn = create_score_function("ref2015_memb")

# Setup FastRelax
relax = FastRelax(scorefxn, 3)

print("Running FastRelax with membrane scoring (3 cycles)...")
relax.apply(pose)

# Get final axis Z-component
ca_a_final = []
ca_b_final = []
for i in range(1, pose.total_residue() + 1):
    res = pose.residue(i)
    if not res.is_polymer():
        continue
    if res.chain() == 1:
        ca_a_final.append([res.atom("CA").xyz().x, res.atom("CA").xyz().y, res.atom("CA").xyz().z])
    elif res.chain() == 2:
        ca_b_final.append([res.atom("CA").xyz().x, res.atom("CA").xyz().y, res.atom("CA").xyz().z])

ca_a_final = np.array(ca_a_final)
ca_b_final = np.array(ca_b_final)
centroid_a_final = ca_a_final.mean(axis=0)
centroid_b_final = ca_b_final.mean(axis=0)
axis_final = centroid_b_final - centroid_a_final
axis_norm_final = axis_final / np.linalg.norm(axis_final)
z_final = abs(axis_norm_final[2])

print(f"Final axis Z-component: {z_final:.4f}")
print(f"Improvement: {z_final - z_initial:.4f}")

# Output relaxed structure
output_pdb = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/main_pipeline/CytbX_4tool/cycle_2/id11_relaxed_membrane.pdb"
pose.dump_pdb(output_pdb)

print(f"Relaxed structure saved to: {output_pdb}")

