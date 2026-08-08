import sys
from Bio.PDB import MMCIFParser, PDBParser
import numpy as np

def get_structure(path):
    parser = MMCIFParser(QUIET=True) if path.endswith('.cif') else PDBParser(QUIET=True)
    return parser.get_structure('s', path)[0]

def burial_ratio(structure, protein_chain, cofactor_chain, cofactor_resn, reference_chain=None, reference_resn=None, radius=8.0):
    protein_atoms = [a.coord for a in structure[protein_chain].get_atoms()]
    cofactor_atoms = [a.coord for a in structure[cofactor_chain].get_atoms() if a.get_parent().resname.strip() == cofactor_resn]

    if not cofactor_atoms:
        return None

    protein_atoms = np.array(protein_atoms)
    cofactor_atoms = np.array(cofactor_atoms)

    count = 0
    for p in protein_atoms:
        dists = np.linalg.norm(cofactor_atoms - p, axis=1)
        if dists.min() <= radius:
            count += 1

    result = {"cofactor_burial_atoms": count, "total_protein_atoms": len(protein_atoms)}

    if reference_chain and reference_resn:
        ref_atoms = np.array([a.coord for a in structure[reference_chain].get_atoms() if a.get_parent().resname.strip() == reference_resn])
        if len(ref_atoms) > 0:
            ref_count = 0
            for p in protein_atoms:
                dists = np.linalg.norm(ref_atoms - p, axis=1)
                if dists.min() <= radius:
                    ref_count += 1
            result["reference_burial_atoms"] = ref_count
            result["burial_ratio_vs_reference"] = round(count / ref_count, 3) if ref_count > 0 else None

    return result

if __name__ == "__main__":
    path = sys.argv[1]
    protein_chain = sys.argv[2] if len(sys.argv) > 2 else "A"
    cofactor_chain = sys.argv[3] if len(sys.argv) > 3 else "C"
    cofactor_resn = sys.argv[4] if len(sys.argv) > 4 else "FMN"
    reference_chain = sys.argv[5] if len(sys.argv) > 5 else "B"
    reference_resn = sys.argv[6] if len(sys.argv) > 6 else "HEM"

    structure = get_structure(path)
    result = burial_ratio(structure, protein_chain, cofactor_chain, cofactor_resn, reference_chain, reference_resn)
    print(f"File: {path}")
    for k, v in result.items():
        print(f"  {k}: {v}")
