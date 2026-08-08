from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley
import itertools
import numpy as np
import sys

def residue_min_distance(res1, res2, ignore_hydrogens=True):
    """Minimum atom-atom distance between two residues."""
    coords1 = []
    coords2 = []
    for atom in res1.get_atoms():
        if ignore_hydrogens and atom.get_name().startswith("H"):
            continue
        coords1.append(atom.coord)
    for atom in res2.get_atoms():
        if ignore_hydrogens and atom.get_name().startswith("H"):
            continue
        coords2.append(atom.coord)
    if not coords1 or not coords2:
        return np.inf
    coords1 = np.asarray(coords1)
    coords2 = np.asarray(coords2)
    dists = np.linalg.norm(coords1[:, None, :] - coords2[None, :, :], axis=2)
    return dists.min()

def chains_in_contact(chain1, chain2, cutoff=8.0):
    """Return True if any residue in chain1 is within cutoff Å of any residue in chain2."""
    residues1 = [r for r in chain1.get_residues() if is_aa(r, standard=True) and r.id[0] == " "]
    residues2 = [r for r in chain2.get_residues() if is_aa(r, standard=True) and r.id[0] == " "]
    if not residues1 or not residues2:
        return False
    for r1 in residues1:
        for r2 in residues2:
            if residue_min_distance(r1, r2) <= cutoff:
                return True
    return False

def compute_chain_sasa(chain, sr):
    """Compute SASA for a single chain (residues only)."""
    sr.compute(chain, level="R")
    return sum(r.sasa for r in chain.get_residues() if is_aa(r, standard=True) and r.id[0] == " ")

def compute_complex_sasa(chains, sr):
    """Compute SASA for a complex of chains."""
    model = chains[0].get_parent()
    sr.compute(model, level="R")
    sasa = 0.0
    for chain in chains:
        sasa += sum(r.sasa for r in chain.get_residues() if is_aa(r, standard=True) and r.id[0] == " ")
    return sasa

def interface_surface_area(structure_file, model_index=0, probe_radius=1.4, n_points=960, contact_cutoff=5.0):
    """Compute pairwise interface SASA between physically contacting chains."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("prot", structure_file)
    model = structure[model_index]
    chains = list(model.get_chains())
    sr = ShrakeRupley(probe_radius=probe_radius, n_points=n_points)

    # Precompute SASA for each chain
    chain_sasa = {chain.id: compute_chain_sasa(chain, sr) for chain in chains}

    results = []

    for chain_a, chain_b in itertools.combinations(chains, 2):
        # Skip if chains are not in contact
        if not chains_in_contact(chain_a, chain_b, cutoff=contact_cutoff):
            continue

        sasa_a = chain_sasa[chain_a.id]
        sasa_b = chain_sasa[chain_b.id]
        sasa_ab = compute_complex_sasa([chain_a, chain_b], sr)

        interface_area = (sasa_a + sasa_b - sasa_ab) / 2.0
        normalised_interface = interface_area / (sasa_a + sasa_b)

        results.append({
            "chain_a": chain_a.id,
            "chain_b": chain_b.id,
            "interface_area": interface_area,
            "normalised_interface": normalised_interface
        })

    return results

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python surface_area_biopython.py <structure.cif>")
        sys.exit(1)

    structure_file = sys.argv[1]
    print(f"Analysing {structure_file}")

    results = interface_surface_area(structure_file)

    outname = structure_file.replace(".cif", "_interface.txt")
    with open(outname, "w") as f:
        f.write("# chainA chainB interface_area normalised_interface\n")
        for r in results:
            f.write(
                f"{r['chain_a']} {r['chain_b']} "
                f"{r['interface_area']:.2f} "
                f"{r['normalised_interface']:.4f}\n"
            )
