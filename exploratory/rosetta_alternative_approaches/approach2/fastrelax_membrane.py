import pyrosetta
import sys

input_pdb = sys.argv[1]
output_pdb = sys.argv[2]
spanfile = sys.argv[3]

pyrosetta.init('-mute all -mp:lipids:has_pore false')

pose = pyrosetta.pose_from_pdb(input_pdb)

from pyrosetta.rosetta.protocols.membrane import AddMembraneMover
amm = AddMembraneMover(spanfile)
amm.apply(pose)

from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
sfxn = ScoreFunctionFactory.create_score_function('franklin2019')

from pyrosetta.rosetta.protocols.relax import FastRelax
fr = FastRelax()
fr.set_scorefxn(sfxn)
fr.constrain_relax_to_start_coords(True)
fr.apply(pose)

pose.dump_pdb(output_pdb)
print(f'Done. Output: {output_pdb}')
print(f'Final score: {sfxn(pose):.3f}')
