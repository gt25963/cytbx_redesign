#!/usr/bin/env python3

# Combine Boltz-2 + ESM3 scores into a single ranked CSV.
# for any cross-state/cross-cycle comparison use rescore_simple.py instead 

import csv
import sys

def normalise(d):
    if not d: return d
    mn, mx = min(d.values()), max(d.values())
    if mx == mn: return {k: 1.0 for k in d}
    return {k: (v - mn) / (mx - mn) for k, v in d.items()}

def main():
    boltz_file = sys.argv[1]
    esm3_file = sys.argv[2]
    output_file = sys.argv[3]

    boltz_scores = {}
    with open(boltz_file) as f:
        for line in f:
            parts = line.strip().split(',')
            if parts[0] == 'id' or len(parts) < 4:
                continue
            try:
                boltz_scores[parts[0]] = float(parts[3])
            except ValueError:
                continue

    esm3_scores = {}
    with open(esm3_file) as f:
        for row in csv.DictReader(f):
            esm3_scores[row["id"]] = float(row["ptm"])

    boltz_norm = normalise(boltz_scores)
    esm3_norm = normalise(esm3_scores)

    all_ids = set(boltz_norm.keys()) | set(esm3_norm.keys())
    combined = {}
    for id in all_ids:
        b = boltz_norm.get(id, 0.0)
        e = esm3_norm.get(id, 0.0)
        combined[id] = (b + e) / 2.0

    sorted_ids = sorted(combined, key=combined.get, reverse=True)

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["id", "boltz_score", "esm3_ptm", "combined_score"])
        for id in sorted_ids:
            writer.writerow([
                id,
                round(boltz_scores.get(id, 0.0), 4),
                round(esm3_scores.get(id, 0.0), 4),
                round(combined[id], 4)
            ])

    print(f"combined scores written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python combine_scores.py <boltz_scores.csv> <esm3_scores.csv> <output.csv>")
        sys.exit(1)
    main()
