import pyrosetta, sys, os
input_pdb  = sys.argv[1]
output_dir = sys.argv[2]
spanfile   = sys.argv[3]
n_decoys   = int(sys.argv[4])
os.makedirs(output_dir, exist_ok=True)

pyrosetta.init('-mute all -mp:lipids:has_pore false '
               '-relax:constrain_relax_to_start_coords '
               '-relax:coord_cst_stdev 1.0 '
               '-extra_res_fa /scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/rq2/ligand_params/FMN.params') ## only needs FMN.params (retained Hem2 uses Rosetta's native HEM residue type)

from pyrosetta.rosetta.protocols.membrane import AddMembraneMover
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    InitializeFromCommandline, OperateOnResidueSubset, RestrictToRepackingRLT)
from pyrosetta.rosetta.core.select.residue_selector import (
    ResidueIndexSelector, NotResidueSelector)
from pyrosetta.rosetta.protocols.relax import FastRelax

# spanfile format: header + antiparallel/n2c + tab-separated spans
pose0 = pyrosetta.pose_from_pdb(input_pdb)
AddMembraneMover(spanfile).apply(pose0)
sfxn = ScoreFunctionFactory.create_score_function('franklin2019')

# design only the HEM_B pocket positions; repack (not design) everything else
design_res = [16,20,23,24,41,74,75,78,81,92,96,99]
design_sel = ResidueIndexSelector(','.join(str(r) for r in design_res))
not_design = NotResidueSelector(design_sel)

print(f"score before: {sfxn(pose0):.1f}")
best = None
for i in range(n_decoys):
    # clone fresh from pose0 each decoy so relax trajectories do NOT compound across decoys
    wp = pose0.clone()
    tf = TaskFactory()
    tf.push_back(InitializeFromCommandline())
    tf.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), not_design))
    fr = FastRelax(sfxn, 3)
    fr.set_task_factory(tf)
    fr.constrain_relax_to_start_coords(True)
    fr.apply(wp)
    sc = sfxn(wp)
    out = os.path.join(output_dir, f'fmnhemb_design_{i:02d}.pdb')
    wp.dump_pdb(out)
    print(f"decoy {i:02d}: score {sc:.1f} -> {out}"); sys.stdout.flush()
    if best is None or sc < best[0]: ## lower REU = better in Rosetta 
        best = (sc, out)
print(f"best: {best[1]} (score {best[0]:.1f})")
