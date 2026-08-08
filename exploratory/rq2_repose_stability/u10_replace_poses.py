#!/usr/bin/env python3
"""
u10_replace_poses.py  (RQ2, U10 track - re-placement)
"""

import argparse
import os
import numpy as np

from fmn_common import (read_atoms, atoms_by_name, centroid,
                        best_fit_plane_normal, apply_transform,
                        write_transformed_ligand)

HEAD_ATOMS = ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9",
              "O1", "O2", "O3", "O4"]
HBOND_EDGE_ATOMS = ["O1", "O2"]


def rotation_between(a, b):
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    v = np.cross(a, b); s = np.linalg.norm(v); c = np.dot(a, b)
    if s < 1e-8:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s ** 2))


def rodrigues(axis, angle_rad):
    k = axis / np.linalg.norm(axis)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(angle_rad) * K + (1 - np.cos(angle_rad)) * (K @ K)


def pose1_geometric(lig_atoms, hemc_coords):
    lig = atoms_by_name(lig_atoms)
    head = [lig[n] for n in HEAD_ATOMS if n in lig]
    head_c = centroid(head)
    head_n = best_fit_plane_normal(head)
    hem_c = centroid(hemc_coords)
    hem_n = best_fit_plane_normal(hemc_coords)
    R = rotation_between(head_n, hem_n)
    t = hem_c - R @ head_c
    return R, t


def sweep_best_pose(lig_atoms, hemc_coords, clash_ca, angle_step_deg=5,
                    include_flip=True):
    lig = atoms_by_name(lig_atoms)
    lig_names = list(lig.keys())
    lig_coords0 = np.array([lig[n] for n in lig_names])

    head = [lig[n] for n in HEAD_ATOMS if n in lig]
    head_c0 = centroid(head)

    R0, t0 = pose1_geometric(lig_atoms, hemc_coords)
    head_c_target = apply_transform([head_c0], R0, t0)[0]
    head_n = best_fit_plane_normal(apply_transform(head, R0, t0))

    arbitrary = np.array([1.0, 0.0, 0.0])
    inplane = arbitrary - np.dot(arbitrary, head_n) * head_n
    if np.linalg.norm(inplane) < 1e-6:
        arbitrary = np.array([0.0, 1.0, 0.0])
        inplane = arbitrary - np.dot(arbitrary, head_n) * head_n
    flip_axis = inplane / np.linalg.norm(inplane)

    def score_for(R, t):
        moved = apply_transform(lig_coords0, R, t)
        worst = float("inf")
        for label, ca_xyz in clash_ca:
            d = min(np.linalg.norm(ca_xyz - m) for m in moved)
            worst = min(worst, d)
        return worst

    flips = [False, True] if include_flip else [False]
    results = []
    best_score, best_R, best_t, best_deg, best_flip = -1.0, None, None, None, None

    for flip in flips:
        Rflip = rodrigues(flip_axis, np.pi) if flip else np.eye(3)
        for deg in range(0, 360, angle_step_deg):
            Rrot = rodrigues(head_n, np.radians(deg))
            R = Rrot @ Rflip @ R0
            t = head_c_target - R @ head_c0
            s = score_for(R, t)
            results.append((deg, flip, s))
            if s > best_score:
                best_score, best_R, best_t, best_deg, best_flip = s, R, t, deg, flip

    return best_R, best_t, best_score, best_deg, best_flip, results


def sweep_translate(lig_atoms, hemc_coords, clash_ca, base_deg, base_flip,
                    extent=1.5, step=0.25):
    lig = atoms_by_name(lig_atoms)
    lig_coords0 = np.array([lig[n] for n in lig.keys()])

    head = [lig[n] for n in HEAD_ATOMS if n in lig]
    head_c0 = centroid(head)

    R0, t0 = pose1_geometric(lig_atoms, hemc_coords)
    head_c_target = apply_transform([head_c0], R0, t0)[0]
    head_n = best_fit_plane_normal(apply_transform(head, R0, t0))

    arbitrary = np.array([1.0, 0.0, 0.0])
    inplane = arbitrary - np.dot(arbitrary, head_n) * head_n
    if np.linalg.norm(inplane) < 1e-6:
        arbitrary = np.array([0.0, 1.0, 0.0])
        inplane = arbitrary - np.dot(arbitrary, head_n) * head_n
    flip_axis = inplane / np.linalg.norm(inplane)

    Rflip = rodrigues(flip_axis, np.pi) if base_flip else np.eye(3)
    Rrot = rodrigues(head_n, np.radians(base_deg))
    R = Rrot @ Rflip @ R0
    base_t = head_c_target - R @ head_c0

    def score_for(t):
        moved = apply_transform(lig_coords0, R, t)
        worst = float("inf")
        for label, ca_xyz in clash_ca:
            d = min(np.linalg.norm(ca_xyz - m) for m in moved)
            worst = min(worst, d)
        return worst

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


def report_clash_distances(out_pdb, clash_chain, clash_residues):
    lig_coords = [a["xyz"] for a in read_atoms(out_pdb, want_resname="U10")]
    ca = {}
    for a in read_atoms(out_pdb, want_chain=clash_chain):
        if a["name"] == "CA":
            ca[a["resnum"]] = a["xyz"]
    print(f"    Ca-to-U10 distances at clash residues:")
    for r in clash_residues:
        if r not in ca:
            print(f"      res{r}: CA not found")
            continue
        d = min(np.linalg.norm(ca[r] - m) for m in lig_coords)
        print(f"      res{r}: {d:.2f} A")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--swap", required=True)
    ap.add_argument("--hemc-ref", required=True)
    ap.add_argument("--hemc-resnum", type=int, default=114)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sweep", action="store_true")
    ap.add_argument("--clash-chain", default="A")
    ap.add_argument("--clash-residues", type=int, nargs="+", default=[9, 67])
    ap.add_argument("--angle-step", type=int, default=5)
    ap.add_argument("--show-all-above", type=float, default=None)
    ap.add_argument("--sweep-translate", action="store_true")
    ap.add_argument("--fixed-angle", type=int, default=None)
    ap.add_argument("--dump-angle", type=int, default=None)
    ap.add_argument("--flip", action="store_true")
    ap.add_argument("--tag", default="custom")
    ap.add_argument("--fixed-flip", action="store_true")
    ap.add_argument("--translate-extent", type=float, default=1.5)
    ap.add_argument("--translate-step", type=float, default=0.25)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    lig_atoms = read_atoms(args.swap, want_resname="U10")
    hemc_coords = [a["xyz"] for a in
                   read_atoms(args.hemc_ref, want_resname="HEM",
                              want_resnum=args.hemc_resnum)]
    if not hemc_coords:
        raise SystemExit(f"No HEM {args.hemc_resnum} found in {args.hemc_ref}")

    print(f"U10 atoms: {len(lig_atoms)} | HEM_C ref atoms: {len(hemc_coords)}")

    if args.dump_angle is not None:
        head = atoms_by_name(lig_atoms)
        head_pts = [head[n] for n in HEAD_ATOMS if n in head]
        head_c0 = centroid(head_pts)
        R0, t0 = pose1_geometric(lig_atoms, hemc_coords)
        head_c_target = apply_transform([head_c0], R0, t0)[0]
        head_n = best_fit_plane_normal(apply_transform(head_pts, R0, t0))

        arbitrary = np.array([1.0, 0.0, 0.0])
        inplane = arbitrary - np.dot(arbitrary, head_n) * head_n
        if np.linalg.norm(inplane) < 1e-6:
            arbitrary = np.array([0.0, 1.0, 0.0])
            inplane = arbitrary - np.dot(arbitrary, head_n) * head_n
        flip_axis = inplane / np.linalg.norm(inplane)

        Rflip = rodrigues(flip_axis, np.pi) if args.flip else np.eye(3)
        Rrot = rodrigues(head_n, np.radians(args.dump_angle))
        R = Rrot @ Rflip @ R0
        t = head_c_target - R @ head_c0

        out = os.path.join(args.outdir, f"u10_pose_{args.tag}.pdb")
        write_transformed_ligand(args.swap, out, "U10", R, t)
        print(f"[dump-angle] angle={args.dump_angle} flip={args.flip} -> {out}")
        report_clash_distances(out, args.clash_chain, args.clash_residues)
        return

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
            lig_atoms, hemc_coords, clash_ca,
            args.fixed_angle, args.fixed_flip,
            extent=args.translate_extent, step=args.translate_step)

        out = os.path.join(args.outdir, "u10_pose_translated.pdb")
        write_transformed_ligand(args.swap, out, "U10", R, best_t)
        print(f"\n[translate-sweep] base angle={args.fixed_angle} flip={args.fixed_flip}")
        print(f"[translate-sweep] best offset dx,dy,dz={off}, "
              f"worst-case Ca-to-U10 distance={score:.2f} A")
        print(f"[translate-sweep] -> {out}")
        report_clash_distances(out, args.clash_chain, args.clash_residues)

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
            lig_atoms, hemc_coords, clash_ca,
            angle_step_deg=args.angle_step, include_flip=True)

        out = os.path.join(args.outdir, "u10_pose_swept.pdb")
        write_transformed_ligand(args.swap, out, "U10", R, t)
        print(f"\n[sweep] best angle={deg} deg, flip={flip}, "
              f"worst-case Ca-to-U10 distance={score:.2f} A")
        print(f"[sweep] -> {out}")
        report_clash_distances(out, args.clash_chain, args.clash_residues)

        if args.show_all_above is not None:
            passing = [r for r in results if r[2] >= args.show_all_above]
            passing.sort(key=lambda r: r[0])
            print(f"\n  all angles clearing {args.show_all_above:.1f} A "
                  f"({len(passing)}/{len(results)} candidates):")
            for deg_r, flip_r, s_r in passing:
                print(f"    angle={deg_r:3d} flip={flip_r}  worst_d={s_r:.2f} A")
        else:
            results.sort(key=lambda r: -r[2])
            print("\n  top 5 candidates (angle, flip, worst-case distance):")
            for deg_r, flip_r, s_r in results[:5]:
                print(f"    angle={deg_r:3d} flip={flip_r}  worst_d={s_r:.2f} A")
        return

    R1, t1 = pose1_geometric(lig_atoms, hemc_coords)
    out1 = os.path.join(args.outdir, "u10_pose1_geometric.pdb")
    write_transformed_ligand(args.swap, out1, "U10", R1, t1)
    print(f"[pose1 geometric] -> {out1}")
    report_clash_distances(out1, args.clash_chain, args.clash_residues)


if __name__ == "__main__":
    main()
