#!/usr/bin/env python3
"""RQ2 aggregate: protein-cofactor interface ipTM.
Chain layout (verified): 0=protein, 1=HEM_B retained, 2=swapped cofactor.
Primary metric = protein-cofactor pair [0][2] symmetrised, max over 5 models.
Also reports protein-HEM [0][1] as a sanity column (native site should stay high).
Usage: python build_rq2_aggregate.py <chai_outputs_dir> <af3_dir> <out_csv>"""
import numpy as np, glob, os, re, csv, json, sys
chai_dir, af3_dir, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]
def bare(s):
    m = re.search(r"id(\d+)", s); return m.group(1) if m else None
def sym(mm, i, j): return float((mm[i,j] + mm[j,i]) / 2)
# --- Chai ---
chai_cof, chai_hem = {}, {}
for d in sorted(glob.glob(f"{chai_dir}/*/")):
    sid = bare(os.path.basename(d.rstrip("/")))
    if not sid: continue
    cof, hem = [], []
    for npz in glob.glob(f"{d}scores.model_idx_*.npz"):
        a = np.load(npz); m = a["per_chain_pair_iptm"]; mm = m[0] if m.ndim==3 else m
        if mm.shape[0] < 3: continue
        cof.append(sym(mm,0,2)); hem.append(sym(mm,0,1))
    if cof:
        chai_cof[sid] = max(cof)
        chai_hem[sid] = max(hem)
# --- AF3: protein-cofactor [0][2] ---
af3_cof = {}
for f in glob.glob(f"{af3_dir}/**/*summary_confidences.json", recursive=True):
    sid = bare(f)
    if not sid: continue
    try:
        j = json.load(open(f)); cm = np.array(j["chain_pair_iptm"])
        if cm.shape[0] < 3: continue
        v = sym(cm,0,2)
        if sid not in af3_cof or v > af3_cof[sid]: af3_cof[sid] = v
    except Exception: pass
ids = sorted(set(chai_cof)|set(af3_cof), key=lambda x:int(x))
rows = []
for sid in ids:
    c = chai_cof.get(sid,0.0); a = af3_cof.get(sid,0.0); h = chai_hem.get(sid,0.0)
    rows.append({"id":sid,"chai_cofactor_iptm":c,"af3_cofactor_iptm":a,
                 "af3_chai_mean":(a+c)/2,"chai_hemB_iptm_sanity":h})
rows.sort(key=lambda r:-r["chai_cofactor_iptm"])
with open(out_csv,"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["id","chai_cofactor_iptm","af3_cofactor_iptm","af3_chai_mean","chai_hemB_iptm_sanity"]); w.writeheader()
    for r in rows: w.writerow({k:(f"{r[k]:.4f}" if k!="id" else r[k]) for k in r})
print(f"wrote {out_csv} ({len(rows)} designs)")
cv=sorted(chai_cof.values()); av=sorted(af3_cof.values()); hv=sorted(chai_hem.values())
if cv: print(f"chai cofactor-iptm: {cv[0]:.4f}..{cv[-1]:.4f}")
if av: print(f"af3 cofactor-iptm:  {av[0]:.4f}..{av[-1]:.4f}")
if hv: print(f"chai HEM_B sanity:  {hv[0]:.4f}..{hv[-1]:.4f}  (should stay high)")
print("\n=== top 10 by Chai protein-cofactor ipTM ===")
print(f"{'id':>6} {'cofac':>7} {'af3':>7} {'mean':>7} {'hemB':>7}")
for r in rows[:10]:
    print(f"{r['id']:>6} {r['chai_cofactor_iptm']:>7.4f} {r['af3_cofactor_iptm']:>7.4f} {r['af3_chai_mean']:>7.4f} {r['chai_hemB_iptm_sanity']:>7.4f}")
