#!/usr/bin/env python3

# aggregate_chai_scores_v2.py - CORRECTED, uses per_chain_pair_iptm from .npz
# given a directory of Chai-1 .npz outputs, computes the true protein-protein ipTM (not the inflated aggregate score) per design and writes a ranked CSV.

import argparse
import re
import sys
from pathlib import Path

import numpy as np


def find_npz_files(root: Path):
    return sorted(root.rglob("scores.model_idx_*.npz"))

def parse_id_from_path(path: Path):
    # design ID - the trailing 'idN' in the parent directory name
    parent = path.parent.name
    m = re.search(r"id(\d+)$", parent)
    if m:
        return m.group(1)
    # fallback - go up path components looking at the layout path
    for part in path.parts[::-1]:
        m = re.search(r"id(\d+)$", part)
        if m:
            return m.group(1)
    # last resort - use raw parent folder name as the id - not ideal but better than it crashing 
    return parent

def protein_pair_iptm(npz_path: Path, protein_chain_indices=(0, 1)):
    d = np.load(npz_path)
    if "per_chain_pair_iptm" not in d:
        return None
    mat = d["per_chain_pair_iptm"]
    if mat.ndim == 3:
        mat = mat[0]
    i, j = protein_chain_indices
    if i >= mat.shape[0] or j >= mat.shape[1]:
        return None
    # average [i][j] and [j][i] since the matrix isn't guaranteed symmetric
    return float((mat[i, j] + mat[j, i]) / 2.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--expected-ids", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=Path("chai_ranking_cycle1_corrected.csv"))
    # default chains (0,1) assumes a 2-chain protein-protein interface (e.g. C2 dimer); override for other states/tracks via --protein-chains
    ap.add_argument("--protein-chains", type=int, nargs=2, default=[0, 1])
    args = ap.parse_args()

    files = find_npz_files(args.root)
    if not files:
        sys.exit(f"No scores.model_idx_*.npz found under {args.root}")

    results = {}
    empty = []
    for f in files:
        sid = parse_id_from_path(f)
        score = protein_pair_iptm(f, tuple(args.protein_chains))
        if score is None:
            empty.append((sid, str(f)))
            continue
        # chai predicts 5 models per design - keep only the best scoring one per id
        if sid not in results or score > results[sid][0]:
            results[sid] = (score, str(f))
    # completeness check - flag ids that never produced a usable score 
    missing = []
    if args.expected_ids and args.expected_ids.exists():
        expected = [l.strip() for l in args.expected_ids.read_text().splitlines() if l.strip()]
        missing = [e for e in expected if e not in results]

    ranked = sorted(results.items(), key=lambda kv: kv[1][0], reverse=True)
    import csv as csv_mod
    with open(args.out, "w", newline="") as fh:
        w = csv_mod.writer(fh)
        w.writerow(["rank", "id", "chai_protein_pair_iptm", "source_file"])
        for rank, (sid, (score, path)) in enumerate(ranked, 1):
            w.writerow([rank, sid, f"{score:.4f}", path])

    # summary - counts here should match expected before trusting the ranking downstream
    print(f"NPZ files found: {len(files)}")
    print(f"Unique sequence ids: {len(results)}")
    print(f"Files with no per_chain_pair_iptm: {len(empty)}")
    if empty:
        for sid, p in empty[:10]:
            print(f"  EMPTY  id={sid}  {p}")
    if args.expected_ids:
        print(f"Expected ids:       {len(expected)}")
        if missing:
            print(f"MISSING ({len(missing)}): {', '.join(missing)}")
        else:
            print("All expected ids present.")
    print(f"Ranking written to: {args.out}")
    if ranked:
        print("\nTop 10 (corrected protein-protein iptm):")
        for rank, (sid, (score, _)) in enumerate(ranked[:10], 1):
            print(f"  {rank}. id{sid}  {score:.4f}")


if __name__ == "__main__":
    main()
