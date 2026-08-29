# stability_worker.py (RQ2, rq2_repose_stability - exploratory)
# NOTE: Part of the earlier, discontinued cofactor stability-scoring approach.
# Measures how far a placed cofactor (FMN or U10) drifts from its starting position after a short Rosetta relax, as a proxy for pose stability.
# Superseded by the LigandMPNN-based redesign pipeline used for the final reported results.

import sys, os, numpy as np, pyrosetta
from pyrosetta import pose_from_file
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory

cofactor = sys.argv[1] ## FMN or U10
design_id = int(sys.argv[2])

work = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline"
params = f"{work}/rq2/ligand_params/{cofactor}.params"
pocket = "FMN_pocket" if cofactor == "FMN" else "U10_pocket"
packed = f"{work}/rq2/design/{pocket}/cycle1_relaxed/LigandMPNN/outputs/packed/holo_hemC_relaxed_v1_packed_{design_id}_1.pdb"
outdir = f"{work}/rq2/design/{pocket}/stability_results"
os.makedirs(outdir, exist_ok=True)
outfile = f"{outdir}/id{design_id}.txt"

# Load the cofactor's params so Rosetta recognises it as a residue type, with sidechain-only relax settings (backbone/coords constrained to start)
pyrosetta.init(f"-mp:lipids:has_pore false -ex1 -ex2 -mute all -extra_res_fa {params} "
               "-relax:constrain_relax_to_start_coords -relax:coord_constrain_sidechains false "
               "-relax:ramp_constraints false")

def cen(p):
    # Centroid of every atom belonging to the cofactor residue, used to measure how far it moves before vs after relax
    c=[]
    for i in range(1,p.total_residue()+1):
        r=p.residue(i)
        if r.name3().strip()==cofactor:
            for a in range(1,r.natoms()+1):
                x=r.xyz(a); c.append([x.x,x.y,x.z])
    return np.mean(np.array(c),axis=0)

# Skip missing input rather than erroring, so a batch of these can run over many design_ids without one missing file killing the whole array job
if not os.path.exists(packed):
    open(outfile,"w").write(f"{design_id},MISSING,NA\n"); sys.exit(0)

try:
    p=pose_from_file(packed)
    sf=ScoreFunctionFactory.create_score_function("ref2015")
    s0=cen(p)
    # 3 rounds of FastRelax, then measure cofactor drift and final energy
    FastRelax(sf,3).apply(p)
    s1=cen(p); e1=sf(p); d=float(np.linalg.norm(s1-s0))
    open(outfile,"w").write(f"{design_id},{d:.3f},{e1:.1f}\n")
except Exception as ex:
    # Catch-all so one bad structure doesn't kill the batch; error is logged per-id instead
    open(outfile,"w").write(f"{design_id},ERROR,{str(ex)[:50]}\n")
