#!/usr/bin/env python3
"""
Interface residue selection for LigandMPNN input.
- Identifies interface residues between subunits within 8A
- Detects haem-coordinating His residues and adds to fixed residues
- Outputs residue list for LigandMPNN --fixed_residues argument
"""

import sys
from Bio.PDB import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch

def get_interface_residues(structure, cutoff=8.0):
    """Find residues at interfaces between different chains."""
    chains = list(structure[0].get_chains())
    all_atoms = list(structure[0].get_atoms())
    ns = NeighborSearch(all_atoms)
    
    interface_residues = set()
    
    for chain_a in chains:
        for chain_b in chains:
            if chain_a.id >= chain_b.id:
                continue
            for residue_a in chain_a.get_residues():
                for atom in residue_a.get_atoms():
                    neighbors = ns.search(atom.coord, cutoff, level='R')
                    for neighbor in neighbors:
                        if neighbor.get_parent().id == chain_b.id:
                            interface_residues.add(
                                (chain_a.id, residue_a.get_id()[1])
                            )
                            interface_residues.add(
                                (chain_b.id, neighbor.get_id()[1])
                            )
    
    return interface_residues


def get_haem_coordinating_residues(structure, ligand_name="HEM", cutoff=6.0):
    """Find His residues coordinating haem iron atoms."""
    coordinating_residues = set()
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip() == ligand_name:
                    for atom in residue:
                        if atom.element == "FE" or atom.name == "FE":
                            all_atoms = list(structure[0].get_atoms())
                            ns = NeighborSearch(all_atoms)
                            neighbors = ns.search(atom.coord, cutoff, level='R')
                            for neighbor in neighbors:
                                if neighbor.get_resname() == "HIS":
                                    chain_id = neighbor.get_parent().id
                                    res_num = neighbor.get_id()[1]
                                    coordinating_residues.add((chain_id, res_num))
    
    return coordinating_residues


def format_for_mpnn(residue_set):
    """Format residue set as LigandMPNN fixed_residues string e.g. A1,A5,B1,B5"""
    sorted_residues = sorted(residue_set, key=lambda x: (x[0], x[1]))
    return " ".join(f"{chain}{resnum}" for chain, resnum in sorted_residues)


def main():
    if len(sys.argv) < 2:
        print("Usage: python biopython_selection.py <input.cif>")
        sys.exit(1)
    
    pdb_path = sys.argv[1]
    

    if pdb_path.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    
    interface_residues = get_interface_residues(structure, cutoff=8.0)
    haem_residues = get_haem_coordinating_residues(structure, ligand_name="HEM", cutoff=6.0)
    
    interface_only = interface_residues - haem_residues
    fixed = haem_residues
    
    print(f"Interface residues (designable): {len(interface_only)}")
    print(f"Haem-coordinating His (fixed): {len(fixed)}")
    
    print(f"INTERFACE={format_for_mpnn(interface_only)}")
    print(f"FIXED={format_for_mpnn(fixed)}")


if __name__ == "__main__":
    main()