#!/usr/bin/env python3
"""
fmn_common.py

Shared helpers for the FMN re-placement / relaxation comparison (RQ2).
Pure-geometry parts use numpy only; Rosetta parts import pyrosetta lazily so
the geometry utilities can be unit-tested without a license check.

Atom-name conventions (from holo_hemC_swap_v2.pdb FMN block):
  Isoalloxazine ring system (the redox-active head, used for superposition):
    benzene ring:      C3 C4 C5 C6 C7   (and C2 bridging)
    pyrazine/pyrimidine: N1 C9 N2 C10 N3 C11 C12 N4
    carbonyls:         O1 (on C10), O2 (on C11)
  H-bonding edge atoms (used for His-anchored pose 2):
    N3 (has H), O1, O2   -- the pyrimidine edge that flavoproteins read
  Tail (NOT used for head placement): C13-C17, hydroxyls O7/O8/O9, phosphate P1/O3-O6
"""

import numpy as np

ISO_RING_ATOMS = ["C2", "C3", "C4", "C5", "C6", "C7",
                  "N1", "C9", "N2", "C10", "N3", "C11", "C12", "N4"]

HBOND_EDGE_ATOMS = ["N3", "O1", "O2"]


def read_atoms(path, want_resname=None, want_resnum=None, want_chain=None,
               heavy_only=True):
    out = []
    for line in open(path):
        if line[:6].strip() not in ("ATOM", "HETATM"):
            continue
        name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21].strip()
        try:
            resnum = int(line[22:26])
        except ValueError:
            continue
        elem = line[76:78].strip()
        if heavy_only and (elem == "H" or (not elem and name.startswith("H"))):
            continue
        if want_resname and resname != want_resname:
            continue
        if want_resnum is not None and resnum != want_resnum:
            continue
        if want_chain and chain != want_chain:
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        out.append({"name": name, "resname": resname, "chain": chain,
                    "resnum": resnum, "xyz": xyz})
    return out


def atoms_by_name(atom_list):
    return {a["name"]: a["xyz"] for a in atom_list}


def centroid(coords):
    return np.mean(np.array(coords), axis=0)


def best_fit_plane_normal(coords):
    pts = np.array(coords)
    c = pts.mean(axis=0)
    _, _, vh = np.linalg.svd(pts - c)
    n = vh[2]
    return n / np.linalg.norm(n)


def kabsch(P, Q):
    P = np.array(P); Q = np.array(Q)
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = Q.mean(axis=0) - R @ P.mean(axis=0)
    return R, t


def apply_transform(coords, R, t):
    return (R @ np.array(coords).T).T + t


def write_transformed_ligand(in_pdb, out_pdb, resname, R, t):
    with open(in_pdb) as fh, open(out_pdb, "w") as out:
        for line in fh:
            if line[:6].strip() in ("ATOM", "HETATM") and line[17:20].strip() == resname:
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                nx, ny, nz = R @ xyz + t
                out.write(f"{line[:30]}{nx:8.3f}{ny:8.3f}{nz:8.3f}{line[54:]}")
            else:
                out.write(line)
