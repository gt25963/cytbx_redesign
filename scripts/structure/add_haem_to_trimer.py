#!/usr/bin/env python3
"""
Add haem ligands from original CytbX monomer to RPXDock trimer output.
Usage: python add_haem_to_trimer.py <rpxdock_trimer.pdb> <cytbx_monomer.pdb> <output.pdb>
"""

import sys
import copy
import numpy as np
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Superimposer import Superimposer


def get_chain_transform(ref_chain, target_chain):
    ref_atoms = [a for a in ref_chain.get_atoms() if a.name == 'CA']
    target_atoms = [a for a in target_chain.get_atoms() if a.name == 'CA']
    min_len = min(len(ref_atoms), len(target_atoms))
    sup = Superimposer()
    sup.set_atoms(target_atoms[:min_len], ref_atoms[:min_len])
    return sup.rotran


def main():
    trimer_path = sys.argv[1]
    monomer_path = sys.argv[2]
    output_path = sys.argv[3]

    parser = PDBParser(QUIET=True)
    trimer = parser.get_structure("trimer", trimer_path)
    monomer = parser.get_structure("monomer", monomer_path)

    # get haem residues from monomer
    haem_residues = [r for r in monomer[0].get_residues()
                     if r.get_resname().strip() == "HEM"]
    print(f"Found {len(haem_residues)} haem ligands in monomer")

    # get monomer protein chain
    mono_chain = list(monomer[0].get_chains())[0]

    # find highest existing residue number to avoid ID clashes
    max_resid = max(r.id[1] for r in trimer[0].get_residues())
    new_haem_id = max_resid + 1

    trimer_chains = list(trimer[0].get_chains())
    print(f"Trimer chains: {[c.id for c in trimer_chains]}")

    for chain in trimer_chains:
        rot, tran = get_chain_transform(mono_chain, chain)
        for haem in haem_residues:
            new_haem = copy.deepcopy(haem)
            # transform coordinates
            for atom in new_haem.get_atoms():
                atom.transform(rot, tran)
            # assign unique residue ID
            new_haem.id = (" ", new_haem_id, " ")
            new_haem.segid = ""
            new_haem_id += 1
            chain.add(new_haem)
            print(f"Added HEM {new_haem_id - 1} to chain {chain.id}")

    io = PDBIO()
    io.set_structure(trimer)
    io.save(output_path)
    print(f"Saved holo trimer to {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python add_haem_to_trimer.py <trimer.pdb> <monomer.pdb> <output.pdb>")
        sys.exit(1)
    main()