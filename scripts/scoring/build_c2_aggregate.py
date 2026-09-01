#!/usr/bin/env python3

#Build highest_aggregate_scores.csv for C2.
# feeds seed selectionm from master_cytbx_4tool

import numpy as np, glob, os, re, csv, json, sys

chai_dir, af3_dir, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]
REDUCE = "max" ## set to max | mean | best_agg

def bare(s):
    m = re.search(r"id(\d+)", s); return m.group(1) if m else None

# Chai: real protein-protein ipTM per design
chai = {}
for d in sorted(glob.glob(f"{chai_dir}/top_scoring.cif_id*/")):
    sid = bare(os.path.basename(d.rstrip("/")))
    if not sid: continue
    per_model = []
    for npz in glob.glob(f"{d}scores.model_idx_*.npz"):
        a = np.load(npz)
        m = a["per_chain_pair_iptm"]
        if m.ndim == 3: m = m[0]
        real_iptm = (m[0,1] + m[1,0]) / 2
        agg = float(a["aggregate_score"][0]) if "aggregate_score" in a else real_iptm
        per_model.append((float(real_iptm), agg))
    if not per_model: continue
    if REDUCE == "max":
        chai[sid] = max(p[0] for p in per_model)
    elif REDUCE == "mean":
        chai[sid] = sum(p[0] for p in per_model) / len(per_model)
    elif REDUCE == "best_agg":
        chai[sid] = max(per_model, key=lambda p: p[1])[0]

# AF3: chain_pair_iptm[0][1] (protein-protein), max over samples
af3 = {}
for f in glob.glob(f"{af3_dir}/**/*summary_confidences.json", recursive=True):
    sid = bare(f)
    if not sid: continue
    try:
        with open(f) as fh: j = json.load(fh)
        v = j["chain_pair_iptm"][0][1]
        if sid not in af3 or v > af3[sid]: af3[sid] = float(v)
    except Exception: pass

ids = sorted(set(chai) | set(af3), key=lambda x: int(x))
rows = []
for sid in ids:
    c = chai.get(sid, 0.0); a = af3.get(sid, 0.0)
    rows.append({"id": sid, "chai_real_iptm": c, "af3_iptm": a,
                 "af3_chai_mean": (a+c)/2})
rows.sort(key=lambda r: -r["chai_real_iptm"])

with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["id","chai_real_iptm","af3_iptm","af3_chai_mean"])
    w.writeheader()
    for r in rows:
        w.writerow({k:(f"{r[k]:.4f}" if k!="id" else r[k]) for k in r})

# summary as sanity check - design counts and score ranges should look reasonable i guess 
print(f"wrote {out_csv} ({len(rows)} designs, REDUCE={REDUCE})")
print(f"chai designs={len(chai)} af3 designs={len(af3)}")
cv = sorted(chai.values()); av = sorted(af3.values())
if cv: print(f"chai real-iptm: {cv[0]:.4f}..{cv[-1]:.4f}")
if av: print(f"af3 iptm:       {av[0]:.4f}..{av[-1]:.4f}")
print("\n=== top 10 by Chai real-ipTM ===")
print(f"{'id':>6} {'chai':>7} {'af3':>7} {'mean':>7}")
for r in rows[:10]:
    print(f"{r['id']:>6} {r['chai_real_iptm']:>7.4f} {r['af3_iptm']:>7.4f} {r['af3_chai_mean']:>7.4f}")
