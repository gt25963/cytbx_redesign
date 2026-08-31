import pyrosetta
from pyrosetta import pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.protocols.relax import FastRelax
pyrosetta.init('-mp:lipids:has_pore false -relax:constrain_relax_to_start_coords '
               '-relax:coord_constrain_sidechains false -relax:coord_cst_stdev 0.5 '
               '-extra_res_fa rq2/ligand_params/FMN.params -mute all')
pose=pose_from_pdb('rq2/design/FMN_hemB_pocket/holo_FMN_hemB_prerelax.pdb')
sf=get_fa_scorefxn()
print('score before:', round(sf(pose),1))
fr=FastRelax(sf, 3) ## standard cycle count
fr.constrain_relax_to_start_coords(True)
fr.apply(pose)
print('score after :', round(sf(pose),1))
pose.dump_pdb('rq2/design/FMN_hemB_pocket/holo_FMN_hemB_relaxed.pdb')
print('wrote holo_FMN_hemB_relaxed.pdb')
