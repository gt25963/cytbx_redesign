#!/usr/bin/env python3

# Build a homo-oligomeric CIF by threading a LigandMPNN sequence onto the original holo PDB backbone.
# Sequence format: chain1:chain2 or chain1:chain2:chain3 (colon-separated, one per chain)

from Bio.PDB import PDBParser, MMCIFIO
from Bio.SeqUtils import seq1
import sys

AA_MAP = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

AA_3LETTER = {v: k for k, v in AA_MAP.items()}

def build_oligomer(holo_pdb, sequence, output_cif):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", holo_pdb)
    
    new_chains = sequence.split(":")
    model = list(structure.get_models())[0]
    protein_chains = [c for c in model.get_chains() 
                      if any(r.id[0] == ' ' for r in c)]
    
    print(f"Structure chains: {len(protein_chains)}, sequence chains: {len(new_chains)}")
    
    for chain, new_seq in zip(protein_chains, new_chains):
        protein_residues = [r for r in chain if r.id[0] == ' ']
        if len(protein_residues) != len(new_seq):
            print(f"ERROR: Chain {chain.id} has {len(protein_residues)} residues but sequence has {len(new_seq)}")
            sys.exit(1)
        for residue, aa in zip(protein_residues, new_seq):
            new_resname = AA_3LETTER.get(aa)
            if new_resname:
                residue.resname = new_resname
        print(f"Chain {chain.id}: threaded {len(new_seq)} residues")
    
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(output_cif)
    print(f"Saved to {output_cif}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python build_oligomer_cif.py <holo.pdb> <sequence> <output.cif>")
        sys.exit(1)
    build_oligomer(sys.argv[1], sys.argv[2], sys.argv[3])
