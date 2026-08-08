#!/usr/bin/env python3
"""
ca_to_ligand_distance.py

RQ2 decision tool. For each flagged residue, measures the distance from its
backbone Calpha to the nearest heavy atom of the substituted cofactor
(FMN or U10), in the original holo_hemC_swap_v2.pdb.

Why: LigandMPNN repacking (pack_with_ligand_context) can change side-chain
identity but cannot move the backbone. If a clash persists after repacking, the
question is whether the Calpha itself sits inside the ligand volume.

Interpretation (heavy-atom-to-Calpha distance d):
    d < 3.0 A   -> BACKBONE problem. Calpha is essentially inside the ligand.
                   No side chain can fix this. Needs backbone relaxation
                   (FastRelax/FastDesign with ligand) or ligand re-placement.
    3.0-4.5 A   -> BORDERLINE. Backbone is close; a small side chain may just
                   fit, but the pocket is cramped. Relaxation likely still
                   needed. Treat as backbone-limited.
    d > 4.5 A   -> ROTAMER-level. Calpha has room; the clash was a side-chain
                   placement issue that repacking should, in principle, address.

The retained haem (resname HEM) is ignored entirely.

Usage:
    python ca_to_ligand_distance.py \
        --pdb design/FMN_pocket/cycle_1/LigandMPNN/inputs/holo_hemC_swap_v2.pdb \
        --ligand FMN --residues 44 67 106 107

    python ca_to_ligand_distance.py \
        --pdb design/U10_pocket/cycle_1/LigandMPNN/inputs/holo_hemC_swap_v2.pdb \
        --ligand U10 --residues 9 67
"""

import argparse
import math


def parse_pdb(path):
    """Return (ca_atoms, ligand_atoms_by_resname).
    ca_atoms: dict (chain, resseq) -> (x,y,z)
    ligand_atoms: dict resname -> list of (x,y,z) for heavy atoms
    """
    ca = {}
    lig = {}
    with open(path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21].strip()
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            element = line[76:78].strip()
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            # skip hydrogens
            if element == "H" or (not element and atom_name.startswith("H")):
                continue
            if rec == "ATOM" and atom_name == "CA":
                ca[(chain, resseq)] = (x, y, z)
            if rec == "HETATM":
                lig.setdefault(resname, []).append((x, y, z))
    return ca, lig


def dist(a, b):
    return math.sqrt(sum((p - q) ** 2 for p, q in zip(a, b)))


def verdict(d):
    if d < 3.0:
        return "BACKBONE  (Calpha inside ligand volume -> relax/re-place)"
    if d <= 4.5:
        return "BORDERLINE (backbone-limited -> relaxation likely needed)"
    return "ROTAMER    (Calpha has room -> side-chain level)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--ligand", required=True, help="ligand resname, e.g. FMN or U10")
    ap.add_argument("--residues", required=True, nargs="+", type=int)
    ap.add_argument("--chain", default=None,
                    help="restrict to this protein chain (default: search all)")
    args = ap.parse_args()

    ca, lig = parse_pdb(args.pdb)
    if args.ligand not in lig:
        raise SystemExit(f"Ligand {args.ligand} not found. Present: {list(lig)}")
    lig_atoms = lig[args.ligand]
    print(f"PDB:    {args.pdb}")
    print(f"Ligand: {args.ligand} ({len(lig_atoms)} heavy atoms)")
    print(f"Other HETATM ignored: "
          f"{', '.join(r for r in lig if r != args.ligand) or 'none'}")
    print()
    print(f"{'res':>5} {'chain':>5} {'min_d(A)':>9}  verdict")
    print("-" * 64)

    backbone_hits = []
    for res in args.residues:
        # find this residue's CA on any chain (or the requested chain)
        candidates = [(ch, rs) for (ch, rs) in ca
                      if rs == res and (args.chain is None or ch == args.chain)]
        if not candidates:
            print(f"{res:>5} {'--':>5} {'NO_CA':>9}  CA not found at this resseq")
            continue
        for (ch, rs) in sorted(candidates):
            ca_xyz = ca[(ch, rs)]
            d = min(dist(ca_xyz, l) for l in lig_atoms)
            v = verdict(d)
            print(f"{rs:>5} {ch:>5} {d:>9.2f}  {v}")
            if d <= 4.5:
                backbone_hits.append((rs, ch, d))

    print("-" * 64)
    if backbone_hits:
        print(f"\n{len(backbone_hits)} position(s) are backbone-limited (d <= 4.5 A).")
        print("=> Backbone relaxation with ligand present, or ligand re-placement,")
        print("   is required before further LigandMPNN cycling is meaningful.")
    else:
        print("\nAll flagged positions are rotamer-level (d > 4.5 A).")
        print("=> Backbone is not the problem; revisit repacking settings.")


if __name__ == "__main__":
    main()
