#!/usr/bin/env python3
"""
Rigid-body re-pose search for an RQ2 cofactor to clear flagged-residue clashes
while keeping the redox ring near its His9/His67 coordination pocket.

Searches rotations (about ring normal + small tilts) x translations applied to
the LIGAND ONLY. Scores worst-case clearance from each flagged residue's CA to
the nearest ring/head-group heavy atom. Keeps the transform that maximizes the
minimum clearance, subject to ring-centroid drift <= max_drift from the original.

Usage:
  python repose_search.py --pdb INPUT.pdb --lig FMN --out OUTPUT.pdb \
      --flagged 44,51,106 --max_drift 3.0
  python repose_search.py --pdb INPUT.pdb --lig U10 --out OUTPUT.pdb \
      --flagged 9,67 --max_drift 3.0
"""
import argparse, itertools, sys
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select

# head-group / redox-ring heavy atoms per cofactor (names from the input PDBs)
RING_ATOMS = {
    "FMN": {"C1","C2","C3","C4","C5","C6","C7","C8","C9","N1","N2","N3","N4"},
    "U10": {"C1","C2","C3","C4","C5","C6","C7","C8","O1","O2","O3","O4"},
}
COORD_HIS = [9, 67]   # His coordinating the swapped-cofactor site

def heavy(atom): return atom.element != "H"

def rotation_matrix(axis, theta):
    axis = axis/np.linalg.norm(axis)
    a = np.cos(theta/2.0)
    b,c,d = -axis*np.sin(theta/2.0)
    return np.array([
        [a*a+b*b-c*c-d*d, 2*(b*c+a*d),     2*(b*d-a*c)],
        [2*(b*c-a*d),     a*a+c*c-b*b-d*d, 2*(c*d+a*b)],
        [2*(b*d+a*c),     2*(c*d-a*b),     a*a+d*d-b*b-c*c]])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--lig", required=True, choices=["FMN","U10"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--flagged", required=True, help="comma-sep residue numbers")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--max_drift", type=float, default=3.0)
    ap.add_argument("--clash_thresh", type=float, default=2.0)
    args = ap.parse_args()

    flagged = [int(x) for x in args.flagged.split(",")]
    s = PDBParser(QUIET=True).get_structure("s", args.pdb)
    model = s[0]
    chain = model[args.chain]

    # collect ligand residue + atoms
    lig_res = None
    for res in chain:
        if res.get_resname().strip() == args.lig:
            lig_res = res; break
    if lig_res is None:
        # ligands sometimes on their own chain
        for ch in model:
            for res in ch:
                if res.get_resname().strip() == args.lig:
                    lig_res = res; break
    if lig_res is None: sys.exit(f"{args.lig} not found")

    lig_atoms = [a for a in lig_res if heavy(a)]
    lig_coords = np.array([a.get_coord() for a in lig_atoms])
    ring_mask = np.array([a.get_name() in RING_ATOMS[args.lig] for a in lig_atoms])
    if ring_mask.sum() < 3: sys.exit("too few ring atoms identified")

    ring0 = lig_coords[ring_mask]
    ring_centroid0 = ring0.mean(axis=0)
    # ring normal via SVD
    u,sv,vt = np.linalg.svd(ring0 - ring_centroid0)
    ring_normal = vt[2]

    # flagged CA coords + coordinating His NE2
    ca = {}
    for r in flagged:
        try: ca[r] = chain[r]["CA"].get_coord()
        except KeyError: pass
    his_ne2 = {}
    for r in COORD_HIS:
        try: his_ne2[r] = chain[r]["NE2"].get_coord()
        except KeyError: pass

    def worst_clearance(coords):
        ring = coords[ring_mask]
        m = np.inf
        for r,cac in ca.items():
            d = np.linalg.norm(ring - cac, axis=1).min()
            m = min(m, d)
        return m

    base_clear = worst_clearance(lig_coords)
    print(f"baseline worst clearance (ring->flagged CA): {base_clear:.3f} A")

    # search grid: rotations about ring normal + two in-plane axes, translations
    in_plane1 = vt[0]; in_plane2 = vt[1]
    angles = np.deg2rad(np.arange(0, 360, 15))
    tilts  = np.deg2rad([-20,-10,0,10,20])
    trans_steps = np.arange(-2.0, 2.01, 1.0)

    best = {"clear": base_clear, "coords": lig_coords.copy(), "desc": "identity",
            "drift": 0.0}

    for a_spin in angles:
        R1 = rotation_matrix(ring_normal, a_spin)
        for a_t1 in tilts:
            R2 = rotation_matrix(in_plane1, a_t1)
            for a_t2 in tilts:
                R3 = rotation_matrix(in_plane2, a_t2)
                R = R3 @ R2 @ R1
                rotated = (lig_coords - ring_centroid0) @ R.T + ring_centroid0
                for tx in trans_steps:
                    for ty in trans_steps:
                        for tz in trans_steps:
                            shift = tx*in_plane1 + ty*in_plane2 + tz*ring_normal
                            cand = rotated + shift
                            new_centroid = cand[ring_mask].mean(axis=0)
                            drift = np.linalg.norm(new_centroid - ring_centroid0)
                            if drift > args.max_drift: continue
                            clr = worst_clearance(cand)
                            if clr > best["clear"]:
                                best = {"clear": clr, "coords": cand.copy(),
                                        "desc": f"spin={np.rad2deg(a_spin):.0f} "
                                                f"tilt=({np.rad2deg(a_t1):.0f},{np.rad2deg(a_t2):.0f}) "
                                                f"trans=({tx},{ty},{tz})",
                                        "drift": drift}
    print(f"best worst-clearance: {best['clear']:.3f} A  (drift {best['drift']:.2f} A)")
    print(f"  transform: {best['desc']}")
    # report His coordination distances at best pose
    ring_best = best["coords"][ring_mask]
    rb_centroid = ring_best.mean(axis=0)
    for r,ne2 in his_ne2.items():
        d = np.linalg.norm(ring_best - ne2, axis=1).min()
        print(f"  His{r}-NE2 -> nearest ring atom: {d:.2f} A")

    # write transformed ligand back
    for atom, newc in zip(lig_atoms, best["coords"]):
        atom.set_coord(newc)
    io = PDBIO(); io.set_structure(s); io.save(args.out)
    print(f"wrote {args.out}")

if __name__ == "__main__":
    main()
