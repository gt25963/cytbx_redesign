#!/usr/bin/env python
"""
Measure FMN displacement between the cycle-1 packed LigandMPNN scaffold
and the AF3 top model for id11, after aligning on protein CA atoms.

Tells us whether AF3 moved the cofactor from the deliberately-constructed
pose. Small displacement => the two candidate seeds are effectively the
same structure (seed from packed, clean provenance). Large displacement
=> AF3 imposed its own geometry (a reason to seed from the controlled
packed pose, not AF3).
"""
import sys
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, Superimposer

packed = "rq2/design/FMN_pocket/cycle1_relaxed/LigandMPNN/outputs/packed/holo_hemC_relaxed_v1_packed_11_1.pdb"
af3    = "rq2/master_pipeline/RQ2_FMN/cycle_1/af3/outputs/batch_2/id11/id11_model.cif"

s_packed = PDBParser(QUIET=True).get_structure("packed", packed)
s_af3    = MMCIFParser(QUIET=True).get_structure("af3", af3)

def ca_by_resseq(structure):
    """Map residue sequence number -> CA atom, protein chain only."""
    d = {}
    for model in structure:
        for chain in model:
            for res in chain:
                if "CA" in res:
                    d[res.id[1]] = res["CA"]
        break
    return d

def fmn_atoms(structure):
    """Return FMN atoms keyed by atom name (first FMN residue found)."""
    out = {}
    for model in structure:
        for chain in model:
            for res in chain:
                if res.resname.strip() == "FMN":
                    for atom in res:
                        out[atom.get_name().strip()] = atom
        break
    return out

ca_p = ca_by_resseq(s_packed)
ca_a = ca_by_resseq(s_af3)
shared = sorted(set(ca_p) & set(ca_a))
print(f"shared CA residues for alignment: {len(shared)}")

fixed  = [ca_p[r] for r in shared]
moving = [ca_a[r] for r in shared]

sup = Superimposer()
sup.set_atoms(fixed, moving)
print(f"protein CA alignment RMSD: {sup.rms:.3f} A")

# apply the protein-derived transform to the whole AF3 structure
sup.apply(list(s_af3.get_atoms()))

fmn_p = fmn_atoms(s_packed)
fmn_a = fmn_atoms(s_af3)
common = sorted(set(fmn_p) & set(fmn_a))
print(f"FMN atoms matched by name: {len(common)}")
if not common:
    print("No FMN atoms matched by name. Atom naming differs; compare by centroid instead.")
    pc = np.mean([a.coord for a in fmn_p.values()], axis=0)
    ac = np.mean([a.coord for a in fmn_a.values()], axis=0)
    print(f"FMN centroid displacement: {np.linalg.norm(pc - ac):.3f} A")
    sys.exit(0)

devs = [np.linalg.norm(fmn_p[n].coord - fmn_a[n].coord) for n in common]
devs = np.array(devs)
print(f"FMN per-atom displacement after protein alignment:")
print(f"  mean {devs.mean():.3f} A, max {devs.max():.3f} A, min {devs.min():.3f} A")

pc = np.mean([fmn_p[n].coord for n in common], axis=0)
ac = np.mean([fmn_a[n].coord for n in common], axis=0)
print(f"FMN centroid displacement: {np.linalg.norm(pc - ac):.3f} A")