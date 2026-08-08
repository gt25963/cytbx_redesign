#!/usr/bin/env python3
"""
Heavy-atom clash check for RQ2 LigandMPNN packed outputs.

Checks whether the originally identified clashing residues still clash
with the ligand in each packed structure, post-LigandMPNN repacking.

FMN track clash residues (chain A): 44, 51, 67, 106, 107
U10 track clash residues (chain A): 9, 67

Usage:
    python check_clashes.py --track fmn --packed_dir /path/to/packed --out results_fmn.csv
    python check_clashes.py --track u10 --packed_dir /path/to/packed --out results_u10.csv

Requires: biopython (pip install biopython --break-system-packages)
"""
import argparse
import csv
import glob
import os
import sys

from Bio.PDB import PDBParser

CLASH_THRESHOLD = 2.0  # Angstroms, heavy atoms only

TRACK_CONFIG = {
    "fmn": {
        "residues": [44, 51, 67, 106, 107],
        "ligand_resname": "FMN",  # explicitly target FMN only -- HEM (retained native haem) excluded
    },
    "u10": {
        "residues": [9, 67],
        "ligand_resname": "U10",  # explicitly target U10 only -- HEM (retained native haem) excluded
    },
    "q8": {
        "residues": [3, 6, 7, 10, 106, 107, 108, 109, 110, 111, 112],
        "ligand_resname": "Q8",  # surface quinone; BOTH hemes (HEM_B and HEM_C) retained and excluded
    },
    "fmn_hemb": {
        "residues": [16, 20, 23, 24, 41, 74, 75, 78, 81, 92, 96, 99],
        "ligand_resname": "FMN",  # FMN at HEM_B site; retained HEM_C excluded
    },
}


def get_heavy_atoms(residue):
    return [atom for atom in residue if atom.element != "H"]


def find_ligand_residues(structure, target_resname):
    """Find residues matching the exact target ligand resname (e.g. FMN or U10).
    Deliberately excludes HEM (the retained native haem) and everything else."""
    ligand_residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                if resname == target_resname:
                    ligand_residues.append(residue)
    return ligand_residues


def check_structure(pdb_path, track, target_residues, chain_id="A"):
    """
    Returns a dict: {resnum: min_distance_to_any_ligand_heavy_atom}
    and overall pass/fail against CLASH_THRESHOLD.
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("s", pdb_path)
    except Exception as e:
        return {"error": str(e)}

    target_resname = TRACK_CONFIG[track]["ligand_resname"]
    ligand_residues = find_ligand_residues(structure, target_resname)
    if not ligand_residues:
        return {"error": f"no_{target_resname}_found"}

    ligand_heavy_atoms = []
    for lig_res in ligand_residues:
        ligand_heavy_atoms.extend(get_heavy_atoms(lig_res))

    if not ligand_heavy_atoms:
        return {"error": "ligand_has_no_heavy_atoms"}

    result = {}
    model = structure[0]
    if chain_id not in model:
        return {"error": f"chain_{chain_id}_not_found"}
    chain = model[chain_id]

    for resnum in target_residues:
        try:
            residue = chain[resnum]
        except KeyError:
            result[resnum] = None  # residue not present (e.g. renumbering issue)
            continue

        min_dist = float("inf")
        closest_lig_atom = None
        closest_res_atom = None
        for res_atom in get_heavy_atoms(residue):
            for lig_atom in ligand_heavy_atoms:
                d = res_atom - lig_atom
                if d < min_dist:
                    min_dist = d
                    closest_lig_atom = lig_atom
                    closest_res_atom = res_atom
        result[resnum] = {
            "min_dist": round(min_dist, 3) if min_dist != float("inf") else None,
            "clash": (min_dist < CLASH_THRESHOLD) if min_dist != float("inf") else None,
            "res_atom": closest_res_atom.get_name() if closest_res_atom else None,
            "lig_atom": closest_lig_atom.get_name() if closest_lig_atom else None,
            "lig_resname": closest_lig_atom.get_parent().get_resname() if closest_lig_atom else None,
        }
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--track", choices=["fmn", "u10", "q8", "fmn_hemb"], required=True)
    ap.add_argument("--packed_dir", required=True, help="Directory containing packed_*.pdb files")
    ap.add_argument("--out", required=True, help="Output CSV path")
    ap.add_argument("--chain", default="A")
    args = ap.parse_args()

    target_residues = TRACK_CONFIG[args.track]["residues"]
    pdb_files = sorted(glob.glob(os.path.join(args.packed_dir, "*.pdb")))

    if not pdb_files:
        print(f"No PDB files found in {args.packed_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Checking {len(pdb_files)} packed structures for track '{args.track}'")
    print(f"Target residues (chain {args.chain}): {target_residues}")
    print(f"Clash threshold: {CLASH_THRESHOLD} A (heavy atoms only)\n")

    fieldnames = ["pdb_file", "overall_clash_free"]
    for resnum in target_residues:
        fieldnames += [f"res{resnum}_min_dist", f"res{resnum}_clash", f"res{resnum}_atoms"]
    fieldnames.append("error")

    n_clean = 0
    n_clash = 0
    n_error = 0

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for pdb_path in pdb_files:
            fname = os.path.basename(pdb_path)
            row = {"pdb_file": fname}
            res = check_structure(pdb_path, args.track, target_residues, chain_id=args.chain)

            if "error" in res:
                row["error"] = res["error"]
                row["overall_clash_free"] = ""
                n_error += 1
                writer.writerow(row)
                continue

            any_clash = False
            any_missing = False
            for resnum in target_residues:
                info = res.get(resnum)
                if info is None:
                    row[f"res{resnum}_min_dist"] = ""
                    row[f"res{resnum}_clash"] = "MISSING"
                    row[f"res{resnum}_atoms"] = ""
                    any_missing = True
                    continue
                row[f"res{resnum}_min_dist"] = info["min_dist"]
                row[f"res{resnum}_clash"] = info["clash"]
                row[f"res{resnum}_atoms"] = f"{info['res_atom']}-{info['lig_resname']}.{info['lig_atom']}"
                if info["clash"]:
                    any_clash = True

            row["overall_clash_free"] = (not any_clash) and (not any_missing)
            row["error"] = "missing_residues" if any_missing else ""

            if row["overall_clash_free"]:
                n_clean += 1
            elif any_clash:
                n_clash += 1

            writer.writerow(row)

    print(f"Done. Results written to {args.out}")
    print(f"  Clash-free structures: {n_clean}")
    print(f"  Structures with clashes: {n_clash}")
    print(f"  Errors/missing data: {n_error}")
    print(f"  Total: {len(pdb_files)}")


if __name__ == "__main__":
    main()
