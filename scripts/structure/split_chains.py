#!/usr/bin/env python3

# Split a single-chain Boltz-2 CIF output into a proper oligomer.

from Bio.PDB import MMCIFParser, MMCIFIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
import sys
import math

def split_chains(input_cif, output_cif, oligomeric_state):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", input_cif)
    
    model = list(structure.get_models())[0]
    
    chain_a = model["A"] ## get all protein residues from chain A
    protein_residues = [r for r in chain_a if r.id[0] == ' ']
    ligand_residues = [r for r in chain_a if r.id[0] != ' ']
    
    total_residues = len(protein_residues)
    chain_length = math.ceil(total_residues / oligomeric_state)
    
    print(f"Total protein residues: {total_residues}")
    print(f"Splitting into {oligomeric_state} chains of ~{chain_length} residues each")
    
    # create new structure
    new_structure = Structure("split")
    new_model = Model(0)
    new_structure.add(new_model)
    
    chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    
    # split protein residues into chains
    for i in range(oligomeric_state):
        new_chain = Chain(chain_ids[i])
        start = i * chain_length
        end = min((i + 1) * chain_length, total_residues)
        residues = protein_residues[start:end]
        
        # renumber residues starting from 1
        for j, residue in enumerate(residues):
            residue.id = (residue.id[0], j + 1, residue.id[2])
            new_chain.add(residue)
        
        new_model.add(new_chain)
        print(f"Chain {chain_ids[i]}: {len(residues)} residues")
    
    # add ligands to a separate chain
    if ligand_residues:
        lig_chain = Chain(chain_ids[oligomeric_state])
        for j, residue in enumerate(ligand_residues):
            residue.id = (residue.id[0], j + 1, residue.id[2])
            lig_chain.add(residue)
        new_model.add(lig_chain)
        print(f"Chain {chain_ids[oligomeric_state]}: {len(ligand_residues)} ligand residues")
    
    io = MMCIFIO()
    io.set_structure(new_structure)
    io.save(output_cif)
    print(f"Saved to {output_cif}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python split_chains.py <input.cif> <output.cif> <oligomeric_state>")
        sys.exit(1)
    split_chains(sys.argv[1], sys.argv[2], int(sys.argv[3]))
