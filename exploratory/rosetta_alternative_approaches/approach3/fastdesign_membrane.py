import pyrosetta
import sys
import os

#Read in the four required arguments: input structure, output directory, spanfile, and how many independent relax decoys to generate
input_pdb = sys.argv[1]
output_dir = sys.argv[2]
spanfile = sys.argv[3]
n_decoys = int(sys.argv[4])
os.makedirs(output_dir, exist_ok=True)

pyrosetta.init('-mute all -mp:lipids:has_pore false')

#Haem-coordinating histidines held fixed throughout, matching the LigandMPNN pipeline's fixed-residue convention
fixed_his = [9, 37, 67, 95, 121, 149, 179, 207]

pose = pyrosetta.pose_from_pdb(input_pdb)

#Apply the membrane environment before scoring/relaxing
from pyrosetta.rosetta.protocols.membrane import AddMembraneMover
amm = AddMembraneMover(spanfile)
amm.apply(pose)

#Membrane-aware score function used throughout Rosetta trials
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
sfxn = ScoreFunctionFactory.create_score_function('franklin2019')

from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    RestrictToRepacking, OperateOnResidueSubset,
    PreventRepackingRLT, RestrictToRepackingRLT
)
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.protocols.relax import FastRelax

#Generate n_decoys independent relaxes from the same starting pose, since FastRelax stochastically explores different local minima each run
scores = []
for i in range(n_decoys):
    work_pose = pose.clone()

    #fix haem-coordinating His, allow repacking elsewhere
    tf = TaskFactory()
    from pyrosetta.rosetta.core.pack.task.operation import InitializeFromCommandline
    tf.push_back(InitializeFromCommandline())
    fixed_selector = ResidueIndexSelector(','.join(str(r) for r in fixed_his))
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), fixed_selector))

    fr = FastRelax()
    fr.set_scorefxn(sfxn)
    fr.set_task_factory(tf)
    fr.apply(work_pose)

    #Score and save this decoy, flushing output so progress is visible live in the SLURM log
    score = sfxn(work_pose)
    out_path = os.path.join(output_dir, f'decoy_{i:04d}.pdb')
    work_pose.dump_pdb(out_path)
    scores.append((score, out_path))
    print(f'Decoy {i:04d}: score {score:.3f}')
    sys.stdout.flush()

#Rank all decoys by score (lower is better) and write a summary file
scores.sort()
with open(os.path.join(output_dir, 'scores.txt'), 'w') as f:
    f.write('rank\tscore\tfile\n')
    for rank, (score, path) in enumerate(scores):
        f.write(f'{rank}\t{score:.3f}\t{os.path.basename(path)}\n')
print(f'Done. Best score: {scores[0][0]:.3f} — {scores[0][1]}')
