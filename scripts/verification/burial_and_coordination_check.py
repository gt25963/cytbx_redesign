import sys, glob, os, csv
from Bio.PDB import MMCIFParser
import numpy as np

parser = MMCIFParser(QUIET=True)

def get_atom_coord(structure, chain, resi, atom_name):
    try:
        for res in structure[chain]:
            if res.id[1] == resi:
                if atom_name in res:
                    return res[atom_name].coord
    except Exception:
        pass
    return None

results = []
for cif in glob.glob("rq2/master_pipeline/RQ2_FMN_hemB/cycle_1/af3/outputs/batch_*/id*/id*_model.cif"):
    design_id = os.path.basename(os.path.dirname(cif)).replace("id", "")
    try:
        structure = parser.get_structure("s", cif)[0]

        protein_atoms = np.array([a.coord for a in structure["A"].get_atoms()])
        fmn_atoms = np.array([a.coord for a in structure["C"].get_atoms() if a.get_parent().resname.strip() == "FMN"])
        hem_atoms = np.array([a.coord for a in structure["B"].get_atoms() if a.get_parent().resname.strip() == "HEM"])

        if len(fmn_atoms) == 0 or len(hem_atoms) == 0:
            continue

        fmn_burial = sum(1 for p in protein_atoms if np.linalg.norm(fmn_atoms - p, axis=1).min() <= 8.0)
        hem_burial = sum(1 for p in protein_atoms if np.linalg.norm(hem_atoms - p, axis=1).min() <= 8.0)
        ratio = fmn_burial / hem_burial if hem_burial > 0 else None

        fe = get_atom_coord(structure, "B", None, "FE")
        fe_coord = None
        for res in structure["B"]:
            if res.resname.strip() == "HEM" and "FE" in res:
                fe_coord = res["FE"].coord
                break

        his37_ne2 = get_atom_coord(structure, "A", 37, "NE2")
        his95_ne2 = get_atom_coord(structure, "A", 95, "NE2")

        his37_dist = np.linalg.norm(his37_ne2 - fe_coord) if his37_ne2 is not None and fe_coord is not None else None
        his95_dist = np.linalg.norm(his95_ne2 - fe_coord) if his95_ne2 is not None and fe_coord is not None else None

        results.append({
            "id": design_id,
            "burial_ratio": round(ratio, 3) if ratio else None,
            "his37_to_fe": round(float(his37_dist), 2) if his37_dist is not None else None,
            "his95_to_fe": round(float(his95_dist), 2) if his95_dist is not None else None,
        })
    except Exception as e:
        print(f"Error on {design_id}: {e}")

results.sort(key=lambda r: -(r["burial_ratio"] or 0))

with open("/tmp/burial_coordination_full.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["id", "burial_ratio", "his37_to_fe", "his95_to_fe"])
    w.writeheader()
    w.writerows(results)

print(f"{'id':>6} {'burial':>8} {'his37_fe':>10} {'his95_fe':>10}")
for r in results:
    print(f"{r['id']:>6} {r['burial_ratio']:>8} {r['his37_to_fe']:>10} {r['his95_to_fe']:>10}")
