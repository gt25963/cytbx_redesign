#!/usr/bin/env python3
# compile_tracks_v2.py
# standalone of the track seed selection - incase i need to reuse or just rerun that step 

import argparse
import csv
from pathlib import Path


def load(path, id_col, score_col):
    d = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            try:
                d[str(r[id_col]).strip()] = float(r[score_col])
            except (ValueError, KeyError, TypeError):
                continue
    return d


def mean(vals):
    vals = [v for v in vals if v is not None]
    return sum(vals) / len(vals) if vals else None


def rank_map(score_dict, ids):
    ordered = sorted(ids, key=lambda s: -score_dict.get(s, float("-inf")))
    return {sid: i + 1 for i, sid in enumerate(ordered)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--af3", required=True, type=Path)
    ap.add_argument("--chai", required=True, type=Path)
    ap.add_argument("--boltz", type=Path, default=None)
    ap.add_argument("--esm", type=Path, default=None)
    ap.add_argument("--id-col", default="id")
    ap.add_argument("--af3-col", default="chain_pair_iptm")
    ap.add_argument("--chai-col", default="chai_protein_pair_iptm")
    ap.add_argument("--boltz-col", default="iptm")
    ap.add_argument("--esm-col", default="ptm")
    ap.add_argument("--out", type=Path, default=Path("tracks_cycle1_v2.csv"))
    args = ap.parse_args()

    af3 = load(args.af3, args.id_col, args.af3_col)
    chai = load(args.chai, args.id_col, args.chai_col)
    boltz = load(args.boltz, args.id_col, args.boltz_col) if args.boltz else {}
    esm = load(args.esm, args.id_col, args.esm_col) if args.esm else {}

    have_boltz = bool(boltz)
    have_esm = bool(esm)

    all_ids = sorted(set(af3) | set(chai), key=lambda x: int(x))
    af3_chai_ids = sorted(set(af3) & set(chai), key=lambda x: int(x))

    af3_rank = rank_map(af3, af3_chai_ids)
    chai_rank = rank_map(chai, af3_chai_ids)

    rows = []
    for sid in all_ids:
        a, c = af3.get(sid), chai.get(sid)
        b = boltz.get(sid)
        e = esm.get(sid)
        t1 = a
        t2 = mean([a, c])
        t3 = mean([a, c, b, e]) if (have_boltz or have_esm) else None
        ar, cr = af3_rank.get(sid), chai_rank.get(sid)
        t4_rank_sum = (ar + cr) if (ar is not None and cr is not None) else None
        t5 = c
        rows.append({
            "id": sid, "af3": a, "chai_corrected": c, "boltz": b, "esm_ptm": e,
            "track1_af3": t1, "track2_af3chai_mean": t2,
            "track3_all4_mean": t3,
            "track4_rank_sum": t4_rank_sum,
            "track5_chai_only": t5,
            "af3_rank": ar, "chai_rank": cr,
        })

    fieldnames = list(rows[0].keys())
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})

    def top(key, reverse=True, n=5):
        scored = [r for r in rows if r[key] is not None]
        return sorted(scored, key=lambda r: r[key], reverse=reverse)[:n]

    print(f"Sequences (AF3 or Chai): {len(rows)}")
    print(f"Sequences in both (used for rank tracks): {len(af3_chai_ids)}")
    if not have_boltz:
        print("NOTE: --boltz not supplied; Track 3 skipped pending the "
              "chain-collision-fixed Boltz-2 rerun.")
    if not have_esm:
        print("NOTE: --esm not supplied; Track 3 partial without it.")
    if have_esm:
        print("NOTE: ESM3 contributes 'ptm' (fold confidence), not an interface-specific metric")

    print("\n TRACK 1 (AF3 only) ")
    for r in top("track1_af3"):
        print(f"  id{r['id']}  {r['track1_af3']:.4f}")

    print("\n TRACK 2 (AF3+Chai mean) ")
    for r in top("track2_af3chai_mean"):
        print(f"  id{r['id']}  {r['track2_af3chai_mean']:.4f}  "
              f"(af3={r['af3']}, chai={r['chai_corrected']})")

    if have_boltz or have_esm:
        print("\n TRACK 3 (all-tool mean)")
        for r in top("track3_all4_mean"):
            print(f"  id{r['id']}  {r['track3_all4_mean']:.4f}")

    print("\n TRACK 4 (rank-sum consensus, lower=better)")
    for r in top("track4_rank_sum", reverse=False):
        print(f"  id{r['id']}  rank_sum={r['track4_rank_sum']}  "
              f"(af3_rank={r['af3_rank']}, chai_rank={r['chai_rank']}, "
              f"af3={r['af3']}, chai={r['chai_corrected']})")

    print("\n TRACK 5 (Chai-1 only, corrected) ")
    for r in top("track5_chai_only"):
        print(f"  id{r['id']}  {r['track5_chai_only']:.4f}  (af3={r['af3']})")

    print(f"\nFull table: {args.out}")


if __name__ == "__main__":
    main()
