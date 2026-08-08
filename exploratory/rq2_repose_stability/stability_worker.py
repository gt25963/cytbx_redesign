import sys, os, numpy as np, pyrosetta
from pyrosetta import pose_from_file
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory

cofactor = sys.argv[1]          # FMN or U10
design_id = int(sys.argv[2])
work = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline"
params = f"{work}/rq2/ligand_params/{cofactor}.params"
pocket = "FMN_pocket" if cofactor == "FMN" else "U10_pocket"
packed = f"{work}/rq2/design/{pocket}/cycle1_relaxed/LigandMPNN/outputs/packed/holo_hemC_relaxed_v1_packed_{design_id}_1.pdb"
outdir = f"{work}/rq2/design/{pocket}/stability_results"
os.makedirs(outdir, exist_ok=True)
outfile = f"{outdir}/id{design_id}.txt"

pyrosetta.init(f"-mp:lipids:has_pore false -ex1 -ex2 -mute all -extra_res_fa {params} "
               "-relax:constrain_relax_to_start_coords -relax:coord_constrain_sidechains false "
               "-relax:ramp_constraints false")

def cen(p):
    c=[]
    for i in range(1,p.total_residue()+1):
        r=p.residue(i)
        if r.name3().strip()==cofactor:
            for a in range(1,r.natoms()+1):
                x=r.xyz(a); c.append([x.x,x.y,x.z])
    return np.mean(np.array(c),axis=0)

if not os.path.exists(packed):
    open(outfile,"w").write(f"{design_id},MISSING,NA\n"); sys.exit(0)
try:
    p=pose_from_file(packed)
    sf=ScoreFunctionFactory.create_score_function("ref2015")
    s0=cen(p)
    FastRelax(sf,3).apply(p)
    s1=cen(p); e1=sf(p); d=float(np.linalg.norm(s1-s0))
    open(outfile,"w").write(f"{design_id},{d:.3f},{e1:.1f}\n")
except Exception as ex:
    open(outfile,"w").write(f"{design_id},ERROR,{str(ex)[:50]}\n")
