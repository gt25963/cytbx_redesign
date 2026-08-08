#!/usr/bin/env python3
"""
aggregate_chai_scores.py
"""

import argparse
import csv
import sys
from pathlib import Path

SCORE_FIELD_CANDIDATES = ["iptm", "iPTM", "iptm_score", "aggregate_score"]


def find_score_files(root: Path):
    return sorted(root.rglob("combined_scores.csv"))


def parse_id_from_path(path: Path):
    import re
    parent = path.parent.name
    m = re.search(r"id(\d+)$", parent)
    if m:
        return m.group(1)
    for part in path.parts[::-1]:
        m = re.search(r"id(\d+)$", part)
        if m:
            return m.group(1)
    return parent


def _parse_bracketed_float(s):
    s = s.strip()
    if s.startswith("[") and s.endswith("]"):
        s = s[1:-1]
    return float(s)


def read_best_score(path: Path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    if not rows:
        return None, None
    headers = rows[0].keys()
    field = next((c for c in SCORE_FIELD_CANDIDATES if c in headers), None)
    if field is None:
        return None, None
    vals = []
    for r in rows:
        try:
            vals.append(_parse_bracketed_float(r[field]))
        except (ValueError, TypeError, KeyError):
            continue
    if not vals:
        return field, None
    return field, max(vals)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--expected-ids", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=Path("chai_ranking_cycle1.csv"))
    args = ap.parse_args()

    files = find_score_files(args.root)
    if not files:
        sys.exit(f"No combined_scores.csv found under {args.root}")

    results = {}
    empty = []
    for f in files:
        sid = parse_id_from_path(f)
        field, score = read_best_score(f)
        if score is None:
            empty.append((sid, str(f)))
            continue
        if sid not in results or score > results[sid][0]:
            results[sid] = (score, field, str(f))

    missing = []
    if args.expected_ids and args.expected_ids.exists():
        expected = [l.strip() for l in args.expected_ids.read_text().splitlines() if l.strip()]
        missing = [e for e in expected if e not in results]

    ranked = sorted(results.items(), key=lambda kv: kv[1][0], reverse=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "id", "chai_iptm", "field_used", "source_file"])
        for rank, (sid, (score, field, path)) in enumerate(ranked, 1):
            w.writerow([rank, sid, f"{score:.4f}", field, path])

    print(f"Files found:        {len(files)}")
    print(f"IDs with a score:   {len(results)}")
    print(f"Files empty/no field: {len(empty)}")
    if empty:
        for sid, p in empty:
            print(f"  EMPTY  id={sid}  {p}")
    if args.expected_ids:
        print(f"Expected ids:       {len(expected)}")
        if missing:
            print(f"MISSING ({len(missing)}): {', '.join(missing)}")
        else:
            print("All expected ids present.")
    print(f"Ranking written to: {args.out}")
    if ranked:
        print("\nTop 5:")
        for rank, (sid, (score, field, _)) in enumerate(ranked[:5], 1):
            print(f"  {rank}. id{sid}  {score:.4f} ({field})")


if __name__ == "__main__":
    main()
