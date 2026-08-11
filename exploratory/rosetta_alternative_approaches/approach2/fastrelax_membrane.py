#Setup
##Import PyRosetta and sys for command-line argument handling
import pyrosetta
import sys

##Read in the three required arguments: input structure, output path, spanfile 
##sys.argv[0] is the script name itself, so real arguments start at index 1
input_pdb = sys.argv[1]
output_pdb = sys.argv[2]
spanfile = sys.argv[3]

##Initialise PyRosetta, muting routine output and disabling the membrane pore-lipid option (not relevant to this single-chain relax)
pyrosetta.init('-mute all -mp:lipids:has_pore false')

##Load the input structure into a PyRosetta pose object
pose = pyrosetta.pose_from_pdb(input_pdb)

#Membrane Setup
##Apply the membrane environment defined by the spanfile, so the score function below correctly accounts for the lipid bilayer
from pyrosetta.rosetta.protocols.membrane import AddMembraneMover
amm = AddMembraneMover(spanfile)
amm.apply(pose)

#Scoring and Relax
##Use the membrane-aware franklin2019 score function, standard for transmembrane protein relax throughout this project
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
sfxn = ScoreFunctionFactory.create_score_function('franklin2019')

##Set up FastRelax, keeping the structure close to its starting coordinates (constrain_relax_to_start_coords) so the relax refines rather than substantially reshapes the input structure
from pyrosetta.rosetta.protocols.relax import FastRelax
fr = FastRelax()
fr.set_scorefxn(sfxn)
fr.constrain_relax_to_start_coords(True)
fr.apply(pose)

#Output
##Save the relaxed structure and report the final score to confirm success
pose.dump_pdb(output_pdb)
print(f'Done. Output: {output_pdb}')
print(f'Final score: {sfxn(pose):.3f}')
