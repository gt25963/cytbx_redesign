#!/usr/bin/env python3
"""
fmn_replace_poses.py  (RQ2, arm B - re-placement)
"""
# NOTE: Part of the earlier, discontinued Rosetta-based FMN-at-Hem1 pocket design attempt. 
# This script tries to geometrically place a swapped-in FMN ligand into the Hem1 site by rigid-body transform (rotation/translation), rather than via LigandMPNN redesign, which is the actual method used for the final reported FMN track. 
# Kept for reference; not part of the pipeline that produced the dissertation's reported results.

import argparse
import os
import numpy as np

from fmn_common import (read_atoms, atoms_by_name, centroid,
                        best_fit_plane_normal, apply_transform,
                        write_transformed_ligand,
                        ISO_RING_ATOMS, HBOND_EDGE_ATOMS)


def rotation_between(a, b):
    # Rotation matrix that maps vector a onto vector b, via Rodrigues' formula
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    v = np.cross(a, b); s = np.linalg.norm(v); c = np.dot(a, b)
    # a and b are parallel/antiparallel: cross product is undefined, handle directly
    if s < 1e-8:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s ** 2))


def rodrigues(axis, angle_rad):
    # Standard Rodrigues' rotation formula: rotate by angle_rad around a given axis
    k = axis / np.linalg.norm(axis)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(angle_rad) * K + (1 - np.cos(angle_rad)) * (K @ K)


def pose1_geometric(fmn_atoms, hemc_coords):
    # Baseline pose: align FMN's isoalloxazine ring plane onto Hem2's ring plane, then translate ring centroid onto haem centroid. 
    # Purely geometric, no consideration of nearby His sidechains yet.
    fmn = atoms_by_name(fmn_atoms)
    ring = [fmn[n] for n in ISO_RING_ATOMS if n in fmn]
    ring_c = centroid(ring)
    ring_n = best_fit_plane_normal(ring)
    hem_c = centroid(hemc_coords)
    hem_n = best_fit_plane_normal(hemc_coords)
    R = rotation_between(ring_n, hem_n)
    t = hem_c - R @ ring_c
    return R, t


def pose2_his_anchored(fmn_atoms, his_n_coords, hemc_coords):
    # Refines pose1 by additionally rotating FMN within its own ring plane so its hydrogen-bonding edge points toward the coordinating His nitrogens, rather than at an arbitrary in-plane orientation
    R0, t0 = pose1_geometric(fmn_atoms, hemc_coords)
    fmn = atoms_by_name(fmn_atoms)
    ring = [fmn[n] for n in ISO_RING_ATOMS if n in fmn]
    edge = [fmn[n] for n in HBOND_EDGE_ATOMS if n in fmn]
    ring_c0 = apply_transform([centroid(ring)], R0, t0)[0]
    edge_c0 = apply_transform([centroid(edge)], R0, t0)[0]
    ring_n = best_fit_plane_normal(apply_transform(ring, R0, t0))
    his_mid = centroid(his_n_coords)

    # Project vectors into the ring plane so the rotation angle is computed purely in-plane (out-of-plane component is irrelevant to this step)
    def in_plane(v):
        return v - np.dot(v, ring_n) * ring_n
    cur = in_plane(edge_c0 - ring_c0)
    tgt = in_plane(his_mid - ring_c0)
    if np.linalg.norm(cur) < 1e-6 or np.linalg.norm(tgt) < 1e-6:
        return R0, t0
    cur /= np.linalg.norm(cur); tgt /= np.linalg.norm(tgt)
    # Signed angle between cur and tgt, using ring_n to determine rotation direction
    cosang = np.clip(np.dot(cur, tgt), -1, 1)
    sinang = np.dot(np.cross(cur, tgt), ring_n)
    ang = np.arctan2(sinang, cosang)
    Rrot = rodrigues(ring_n, ang)
    R = Rrot @ R0
    t = Rrot @ (t0 - ring_c0) + ring_c0
    return R, t

def sweep_best_pose(fmn_atoms, hemc_coords, clash_ca, angle_step_deg=10,
                    include_flip=True):
    # Brute-force search over in-plane rotation angle (and optionally a 180-degree ring flip) to find the orientation that maximises minimum distance to nearby clash-prone Ca atoms - i.e. the least clashing pose
    fmn = atoms_by_name(fmn_atoms)
    fmn_names = list(fmn.keys())
    fmn_coords0 = np.array([fmn[n] for n in fmn_names])

    ring = [fmn[n] for n in ISO_RING_ATOMS if n in fmn]
    ring_c0 = centroid(ring)

    R0, t0 = pose1_geometric(fmn_atoms, hemc_coords)
    ring_c_target = apply_transform([ring_c0], R0, t0)[0]
    ring_n = best_fit_plane_normal(apply_transform(ring, R0, t0))

    # Build an arbitrary axis lying within the ring plane to flip around (falls back to a second candidate if the first happens to be parallel to ring_n)
    arbitrary = np.array([1.0, 0.0, 0.0])
    inplane = arbitrary - np.dot(arbitrary, ring_n) * ring_n
    if np.linalg.norm(inplane) < 1e-6:
        arbitrary = np.array([0.0, 1.0, 0.0])
        inplane = arbitrary - np.dot(arbitrary, ring_n) * ring_n
    flip_axis = inplane / np.linalg.norm(inplane)

    # Score = worst-case (minimum) distance from any FMN atom to any clash Ca - higher is better (further from clashing residues)
    def score_for(R, t):
        moved = apply_transform(fmn_coords0, R, t)
        worst = float("inf")
        for label, ca_xyz in clash_ca:
            d = min(np.linalg.norm(ca_xyz - m) for m in moved)
            worst = min(worst, d)
        return worst

    flips = [False, True] if include_flip else [False]
    results = []
    best_score, best_R, best_t, best_deg, best_flip = -1.0, None, None, None, None

    # Grid search over flip x rotation angle, keeping the single best-scoring pose
    for flip in flips:
        Rflip = rodrigues(flip_axis, np.pi) if flip else np.eye(3)
        for deg in range(0, 360, angle_step_deg):
            Rrot = rodrigues(ring_n, np.radians(deg))
            R = Rrot @ Rflip @ R0
            t = ring_c_target - R @ ring_c0
            s = score_for(R, t)
            results.append((deg, flip, s))
            if s > best_score:
                best_score, best_R, best_t, best_deg, best_flip = s, R, t, deg, flip

    return best_R, best_t, best_score, best_deg, best_flip, results


def sweep_translate(fmn_atoms, hemc_coords, clash_ca, base_R, base_deg,
                    base_flip, extent=1.5, step=0.25):
    # Given a fixed rotation (typically the best one from sweep_best_pose), brute-force search a small 3D translation grid around it for further clash reduction. 
    # base_R argument is unused - rotation is recomputed from base_deg/base_flip instead
    fmn = atoms_by_name(fmn_atoms)
    fmn_coords0 = np.array([fmn[n] for n in fmn.keys()])

    ring = [fmn[n] for n in ISO_RING_ATOMS if n in fmn]
    ring_c0 = centroid(ring)

    R0, t0 = pose1_geometric(fmn_atoms, hemc_coords)
    ring_c_target = apply_transform([ring_c0], R0, t0)[0]
    ring_n = best_fit_plane_normal(apply_transform(ring, R0, t0))

    arbitrary = np.array([1.0, 0.0, 0.0])
    inplane = arbitrary - np.dot(arbitrary, ring_n) * ring_n
    if np.linalg.norm(inplane) < 1e-6:
        arbitrary = np.array([0.0, 1.0, 0.0])
        inplane = arbitrary - np.dot(arbitrary, ring_n) * ring_n
    flip_axis = inplane / np.linalg.norm(inplane)

    # Rebuild the fixed rotation from the chosen angle/flip
    Rflip = rodrigues(flip_axis, np.pi) if base_flip else np.eye(3)
    Rrot = rodrigues(ring_n, np.radians(base_deg))
    R = Rrot @ Rflip @ R0
    base_t = ring_c_target - R @ ring_c0

    def score_for(t):
        moved = apply_transform(fmn_coords0, R, t)
        worst = float("inf")
        for label, ca_xyz in clash_ca:
            d = min(np.linalg.norm(ca_xyz - m) for m in moved)
            worst = min(worst, d)
        return worst

    # Exhaustive dx/dy/dz grid search within +/- extent at the given step size
    offsets = np.arange(-extent, extent + 1e-9, step)
    results = []
    best_score, best_t, best_off = -1.0, None, None
    for dx in offsets:
        for dy in offsets:
            for dz in offsets:
                t = base_t + np.array([dx, dy, dz])
                s = score_for(t)
                results.append((dx, dy, dz, s))
                if s > best_score:
                    best_score, best_t, best_off = s, t, (dx, dy, dz)

    return R, best_t, best_score, best_off, results


def report_edge_to_his(out_pdb, his_n_coords):
    # Prints the FMN hydrogen-bonding edge's distance to each candidate His nitrogen, sorted closest-first, as a quick check of H-bond geometry
    fmn = atoms_by_name(read_atoms(out_pdb, want_resname="FMN"))
    edge = [fmn[n] for n in HBOND_EDGE_ATOMS if n in fmn]
    ec = centroid(edge)
    ds = [np.linalg.norm(ec - h) for h in his_n_coords]
    print(f"    edge-centroid to His N atoms: "
          f"{', '.join(f'{d:.2f}' for d in sorted(ds))} A")


def report_clash_distances(out_pdb, clash_chain, clash_residues):
    # Prints minimum Ca-to-FMN distance at each specified clash-prone residue
    fmn_coords = [a["xyz"] for a in read_atoms(out_pdb, want_resname="FMN")]
    ca = {}
    for a in read_atoms(out_pdb, want_chain=clash_chain):
        if a["name"] == "CA":
            ca[a["resnum"]] = a["xyz"]
    print(f"    Ca-to-FMN distances at clash residues:")
    for r in clash_residues:
        if r not in ca:
            print(f"      res{r}: CA not found")
            continue
        d = min(np.linalg.norm(ca[r] - m) for m in fmn_coords)
        print(f"      res{r}: {d:.2f} A")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--swap", required=True)
    ap.add_argument("--hemc-ref", required=True)
    ap.add_argument("--hemc-resnum", type=int, default=114)
    ap.add_argument("--his", type=int, nargs="+", default=[9, 67])
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sweep", action="store_true")
    ap.add_argument("--clash-chain", default="A")
    ap.add_argument("--clash-residues", type=int, nargs="+",
                    default=[44, 67, 106, 107])
    ap.add_argument("--angle-step", type=int, default=10)
    ap.add_argument("--sweep-translate", action="store_true")
    ap.add_argument("--fixed-angle", type=int, default=None)
    ap.add_argument("--fixed-flip", action="store_true")
    ap.add_argument("--translate-extent", type=float, default=1.5)
    ap.add_argument("--translate-step", type=float, default=0.25)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    fmn_atoms = read_atoms(args.swap, want_resname="FMN")
    hemc_coords = [a["xyz"] for a in
                   read_atoms(args.hemc_ref, want_resname="HEM",
                              want_resnum=args.hemc_resnum)]
    if not hemc_coords:
        raise SystemExit(f"No HEM {args.hemc_resnum} found in {args.hemc_ref}")

    # Collect His ND1/NE2 atoms across all specified coordinating residues
    his_n = []
    for hr in args.his:
        for a in read_atoms(args.swap, want_resnum=hr):
            if a["name"] in ("ND1", "NE2"):
                his_n.append(a["xyz"])
    if not his_n:
        raise SystemExit("No His ND1/NE2 atoms found")

    print(f"FMN atoms: {len(fmn_atoms)} | HEM_C ref atoms: {len(hemc_coords)} | "
          f"His N atoms: {len(his_n)}")

    # Three mutually exclusive modes, in order of precedence: fine-grained translation sweep (needs a fixed angle from a prior --sweep run), rotation/flip sweep, or the default two-pose comparison
    if args.sweep_translate:
        if args.fixed_angle is None:
            raise SystemExit("--sweep-translate requires --fixed-angle")
        clash_ca = []
        for r in args.clash_residues:
            hits = [a for a in read_atoms(args.swap, want_chain=args.clash_chain,
                                          want_resnum=r) if a["name"] == "CA"]
            if not hits:
                print(f"  WARNING: no CA found for res{r} chain {args.clash_chain}")
                continue
            clash_ca.append((f"res{r}", hits[0]["xyz"]))

        R, best_t, score, off, results = sweep_translate(
            fmn_atoms, hemc_coords, clash_ca, None,
            args.fixed_angle, args.fixed_flip,
            extent=args.translate_extent, step=args.translate_step)

        out = os.path.join(args.outdir, "fmn_pose4_translated.pdb")
        write_transformed_ligand(args.swap, out, "FMN", R, best_t)
        print(f"\n[translate-sweep] base angle={args.fixed_angle} flip={args.fixed_flip}")
        print(f"[translate-sweep] best offset dx,dy,dz={off}, "
              f"worst-case Ca-to-FMN distance={score:.2f} A")
        print(f"[translate-sweep] -> {out}")
        report_clash_distances(out, args.clash_chain, args.clash_residues)
        report_edge_to_his(out, his_n)

        results.sort(key=lambda r: -r[3])
        print("\n  top 5 offsets (dx, dy, dz, worst-case distance):")
        for dx, dy, dz, s in results[:5]:
            print(f"    ({dx:+.2f}, {dy:+.2f}, {dz:+.2f})  worst_d={s:.2f} A")
        return

    if args.sweep:
        clash_ca = []
        for r in args.clash_residues:
            hits = [a for a in read_atoms(args.swap, want_chain=args.clash_chain,
                                          want_resnum=r) if a["name"] == "CA"]
            if not hits:
                print(f"  WARNING: no CA found for res{r} chain {args.clash_chain}")
                continue
            clash_ca.append((f"res{r}", hits[0]["xyz"]))

        R, t, score, deg, flip, results = sweep_best_pose(
            fmn_atoms, hemc_coords, clash_ca,
            angle_step_deg=args.angle_step, include_flip=True)

        out = os.path.join(args.outdir, "fmn_pose3_swept.pdb")
        write_transformed_ligand(args.swap, out, "FMN", R, t)
        print(f"\n[sweep] best angle={deg} deg, flip={flip}, "
              f"worst-case Ca-to-FMN distance={score:.2f} A")
        print(f"[sweep] -> {out}")
        report_clash_distances(out, args.clash_chain, args.clash_residues)
        report_edge_to_his(out, his_n)

        results.sort(key=lambda r: -r[2])
        print("\n  top 5 candidates (angle, flip, worst-case distance):")
        for deg_r, flip_r, s_r in results[:5]:
            print(f"    angle={deg_r:3d} flip={flip_r}  worst_d={s_r:.2f} A")
        return

    # Default mode: just write out both baseline poses for direct comparison
    R1, t1 = pose1_geometric(fmn_atoms, hemc_coords)
    out1 = os.path.join(args.outdir, "fmn_pose1_geometric.pdb")
    write_transformed_ligand(args.swap, out1, "FMN", R1, t1)
    print(f"[pose1 geometric] -> {out1}")
    report_edge_to_his(out1, his_n)

    R2, t2 = pose2_his_anchored(fmn_atoms, his_n, hemc_coords)
    out2 = os.path.join(args.outdir, "fmn_pose2_hisanchored.pdb")
    write_transformed_ligand(args.swap, out2, "FMN", R2, t2)
    print(f"[pose2 his-anchored] -> {out2}")
    report_edge_to_his(out2, his_n)


if __name__ == "__main__":
    main()
