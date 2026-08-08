import json, glob, re, csv, numpy as np, sys
base, out = sys.argv[1], sys.argv[2]
def cof(pc):
    return float((pc['0']['2'] + pc['2']['0']) / 2)
rows = {}
for d in glob.glob(f"{base}/*_id*"):
    m = re.search(r"id(\d+)$", d)
    if not m: continue
    sid = m.group(1); best = None
    for jf in glob.glob(f"{d}/confidence_*_model_*.json"):
        pc = json.load(open(jf)).get("pair_chains_iptm")
        if not pc or '2' not in pc: continue
        v = cof(pc)
        if best is None or v > best: best = v
    if best is not None: rows[sid] = best
with open(out,"w",newline="") as f:
    w=csv.writer(f); w.writerow(["id","iptm"])
    for sid in sorted(rows,key=int): w.writerow([sid,f"{rows[sid]:.4f}"])
print(f"wrote {out} ({len(rows)} designs)")
if rows:
    vv=sorted(rows.values()); print(f"boltz cofactor-iptm: {vv[0]:.4f}..{vv[-1]:.4f}")
    print("id118 =", rows.get("118"), " id104 =", rows.get("104"))
else:
    print("STILL 0 - check dir names")
