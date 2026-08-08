import json, glob, re, csv, numpy as np
base = "main_pipeline/CytbX_4tool_C3/cycle_1/boltz/outputs/boltz_results_input/predictions"
N_PROT = 3
def pair3(pc):
    # pc is dict-of-dicts keyed by str chain idx
    vals = []
    for i in range(N_PROT):
        for j in range(i+1, N_PROT):
            a = pc[str(i)][str(j)]; b = pc[str(j)][str(i)]
            vals.append((a+b)/2)
    return float(np.mean(vals))
rows = {}
for d in glob.glob(f"{base}/top_scoring.cif_id*"):
    m = re.search(r"id(\d+)$", d)
    if not m: continue
    sid = m.group(1)
    best = None
    for jf in glob.glob(f"{d}/confidence_*_model_*.json"):
        j = json.load(open(jf))
        pc = j.get("pair_chains_iptm")
        if not pc: continue
        v = pair3(pc)
        if best is None or v > best: best = v
    if best is not None: rows[sid] = best
with open("/tmp/c3_boltz.csv","w",newline="") as f:
    w = csv.writer(f); w.writerow(["id","iptm"])
    for sid in sorted(rows, key=int): w.writerow([sid, f"{rows[sid]:.4f}"])
print(f"wrote /tmp/c3_boltz.csv ({len(rows)} designs)")
vv = sorted(rows.values()); print(f"boltz 3-pair iptm: {vv[0]:.4f}..{vv[-1]:.4f}")
print("id14 =", rows.get("14"), " id107 =", rows.get("107"))
