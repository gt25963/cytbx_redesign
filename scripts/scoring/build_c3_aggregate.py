#!/usr/bin/env python3
"""C3 trimer interface ipTM = mean over ALL protein-protein pairs, symmetrised,
over the first N_PROT chains. Max over 5 Chai models."""
import numpy as np, glob, os, re, csv, json, sys
chai_dir, af3_dir, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]
N_PROT = int(sys.argv[4]) if len(sys.argv) > 4 else 3
def bare(s):
    m = re.search(r"id(\d+)", s); return m.group(1) if m else None
def trimer(mm):
    vals = [(mm[i,j]+mm[j,i])/2 for i in range(N_PROT) for j in range(i+1,N_PROT)]
    return float(np.mean(vals))
chai = {}
for d in sorted(glob.glob(f"{chai_dir}/top_scoring.cif_id*/")):
    sid = bare(os.path.basename(d.rstrip("/")))
    if not sid: continue
    pm = []
    for npz in glob.glob(f"{d}scores.model_idx_*.npz"):
        a = np.load(npz); m = a["per_chain_pair_iptm"]; mm = m[0] if m.ndim==3 else m
        pm.append(trimer(mm))
    if pm: chai[sid] = max(pm)
af3 = {}
for f in glob.glob(f"{af3_dir}/**/*summary_confidences.json", recursive=True):
    sid = bare(f)
    if not sid: continue
    try:
        j = json.load(open(f)); cm = np.array(j["chain_pair_iptm"]); v = trimer(cm)
        if sid not in af3 or v > af3[sid]: af3[sid] = float(v)
    except Exception: pass
ids = sorted(set(chai)|set(af3), key=lambda x:int(x))
rows = [{"id":s,"chai_trimer_iptm":chai.get(s,0.0),"af3_trimer_iptm":af3.get(s,0.0),
         "af3_chai_mean":(chai.get(s,0.0)+af3.get(s,0.0))/2} for s in ids]
rows.sort(key=lambda r:-r["chai_trimer_iptm"])
with open(out_csv,"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["id","chai_trimer_iptm","af3_trimer_iptm","af3_chai_mean"]); w.writeheader()
    for r in rows: w.writerow({k:(f"{r[k]:.4f}" if k!="id" else r[k]) for k in r})
print(f"wrote {out_csv} ({len(rows)} designs, N_PROT={N_PROT})")
cv=sorted(chai.values()); av=sorted(af3.values())
if cv: print(f"chai trimer-iptm: {cv[0]:.4f}..{cv[-1]:.4f}")
if av: print(f"af3 trimer-iptm:  {av[0]:.4f}..{av[-1]:.4f}")
print("\n=== top 10 by Chai trimer-ipTM (all 3 pairs averaged) ===")
print(f"{'id':>6} {'chai':>7} {'af3':>7} {'mean':>7}")
for r in rows[:10]: print(f"{r['id']:>6} {r['chai_trimer_iptm']:>7.4f} {r['af3_trimer_iptm']:>7.4f} {r['af3_chai_mean']:>7.4f}")
