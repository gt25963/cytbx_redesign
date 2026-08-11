#!/usr/bin/env python3
"""
FastRelax id11 with membrane scoring to improve C2 axis alignment.
First convert CIF to PDB and remove ligands.
"""
from pyrosetta import *
from pyrosetta.rosetta.protocols.relax import FastRelax
import numpy as np
import os

#Convert CIF to PDB first (remove ligands), since PyRosetta expects a clean protein-only PDB rather than the raw Chai-1 CIF output with cofactors attached
from Bio.PDB import PDBParser, PDBIO, Select

class ProteinSelect(Select):
    #Only keep standard polymer residues from chains A and B (the c2 dimer's two protein chains), filtering out any ligand/cofactor entries
    def accept_residue(self, residue):
        if residue.get_parent().get_id() in ['A', 'B']:
            return 1
        return 0

cif_file = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/main_pipeline/CytbX_4tool/cycle_2/chai/outputs/top_scoring.cif_id11/pred.model_idx_0.cif"
pdb_file = "/tmp/id11_protein_only.pdb"

#Load the CIF with Bio.PDB and write out a protein-only PDB
from Bio.PDB import MMCIFParser
parser = MMCIFParser()
structure = parser.get_structure("id11", cif_file)
io = PDBIO()
io.set_structure(structure)
io.save(pdb_file, ProteinSelect())
print(f"Converted CIF to PDB (protein only): {pdb_file}")

#Load the cleaned structure into PyRosetta
init(extra_options="-mp:lipids:has_pore false")
pose = pose_from_file(pdb_file)
print(">>> FastRelax id11 with membrane scoring")
print(f"Input pose: {pose.total_residue()} residues")

#Measure the initial C2 symmetry axis alignment with the membrane normal, using each chain's Ca centroid to define the inter-chain axis; a Z-component close to 1 means the axis is well-aligned with the membrane normal
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

#Membrane-aware score function (ref2015_memb, note: distinct from franklin2019 used elsewhere in this project's Rosetta trials)
scorefxn = create_score_function("ref2015_memb")

#Run FastRelax to see whether relaxation itself improves axis alignment
relax = FastRelax(scorefxn, 3)
print("Running FastRelax with membrane scoring (3 cycles)...")
relax.apply(pose)

#Re-measure the axis alignment after relax, using the same centroid method
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

#Save the relaxed structure regardless of whether alignment improved
output_pdb = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/main_pipeline/CytbX_4tool/cycle_2/id11_relaxed_membrane.pdb"
pose.dump_pdb(output_pdb)
print(f"Relaxed structure saved to: {output_pdb}")
