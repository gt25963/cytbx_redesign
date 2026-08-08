#!/usr/bin/env python3
"""C2 aggregate with Chai real protein-protein ipTM (max over 5 models)
and AF3 reduced 3 ways for comparison: top-level aggregate, mean of samples, max of samples.
Usage: python build_c2_aggregate_cmp.py <chai_outputs_dir> <af3_dir> <out_csv>"""
import numpy as np, glob, os, re, csv, json, sys

chai_dir, af3_dir, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]

def bare(s):
    m = re.search(r"id(\d+)", s); return m.group(1) if m else None

# --- Chai: real protein-protein ipTM = avg([0][1],[1][0]), max over 5 models ---
chai = {}
for d in sorted(glob.glob(f"{chai_dir}/top_scoring.cif_id*/")):
    sid = bare(os.path.basename(d.rstrip("/")))
    if not sid: continue
    vals = []
    for npz in glob.glob(f"{d}scores.model_idx_*.npz"):
        m = np.load(npz)["per_chain_pair_iptm"]
        if m.ndim == 3: m = m[0]
        vals.append((m[0,1] + m[1,0]) / 2)
    if vals: chai[sid] = float(max(vals))

# --- AF3: three reductions on chain_pair_iptm[0][1] ---
af3_agg, af3_samples = {}, {}
for f in glob.glob(f"{af3_dir}/batch_*/id*/id*_summary_confidences.json"):
    parent = os.path.basename(os.path.dirname(f))
    sid = bare(os.path.basename(f))
    if not sid: continue
    try:
        v = float(json.load(open(f))["chain_pair_iptm"][0][1])
    except Exception:
        continue
    if "seed-" in parent:                       # individual sample
        af3_samples.setdefault(sid, []).append(v)
    else:                                        # top-level aggregate
        af3_agg[sid] = v

af3_mean = {k: sum(v)/len(v) for k, v in af3_samples.items()}
af3_max  = {k: max(v)        for k, v in af3_samples.items()}

ids = sorted(set(chai) | set(af3_agg) | set(af3_mean), key=lambda x: int(x))
rows = []
for sid in ids:
    c = chai.get(sid, 0.0)
    rows.append({"id": sid, "chai_real_iptm": c,
                 "af3_agg": af3_agg.get(sid, 0.0),
                 "af3_mean": af3_mean.get(sid, 0.0),
                 "af3_max": af3_max.get(sid, 0.0),
                 "mean_chai_af3agg":  (c + af3_agg.get(sid,0.0))/2,
                 "mean_chai_af3mean": (c + af3_mean.get(sid,0.0))/2})
rows.sort(key=lambda r: -r["chai_real_iptm"])

fn = ["id","chai_real_iptm","af3_agg","af3_mean","af3_max","mean_chai_af3agg","mean_chai_af3mean"]
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fn); w.writeheader()
    for r in rows:
        w.writerow({k:(f"{r[k]:.4f}" if k!="id" else r[k]) for k in fn})

def spread(d):
    v = sorted(d.values()); return f"{v[0]:.3f}..{v[-1]:.3f} (range {v[-1]-v[0]:.3f}, sd {np.std(list(d.values())):.3f})" if v else "EMPTY"
print(f"wrote {out_csv} ({len(rows)} designs)")
print(f"chai      : {spread(chai)}")
print(f"af3_agg   : {spread(af3_agg)}")
print(f"af3_mean  : {spread(af3_mean)}")
print(f"af3_max   : {spread(af3_max)}")
print("\n=== top 12 by Chai real-ipTM ===")
print(f"{'id':>5} {'chai':>7} {'af3agg':>7} {'af3mn':>7} {'af3mx':>7}")
for r in rows[:12]:
    print(f"{r['id']:>5} {r['chai_real_iptm']:>7.3f} {r['af3_agg']:>7.3f} {r['af3_mean']:>7.3f} {r['af3_max']:>7.3f}")

# rank correlation: does AF3 reduction agree with Chai ordering?
from math import isnan
def spearman(a, b):
    common = [k for k in a if k in b]
    ra = {k:i for i,k in enumerate(sorted(common, key=lambda k:-a[k]))}
    rb = {k:i for i,k in enumerate(sorted(common, key=lambda k:-b[k]))}
    n = len(common)
    if n < 2: return float('nan')
    dd = sum((ra[k]-rb[k])**2 for k in common)
    return 1 - 6*dd/(n*(n*n-1))
print(f"\nSpearman vs Chai ordering:")
print(f"  af3_agg  : {spearman(chai, af3_agg):+.3f}")
print(f"  af3_mean : {spearman(chai, af3_mean):+.3f}")
print(f"  af3_max  : {spearman(chai, af3_max):+.3f}")
