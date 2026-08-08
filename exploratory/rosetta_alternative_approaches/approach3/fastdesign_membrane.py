import pyrosetta
import sys
import os

input_pdb = sys.argv[1]
output_dir = sys.argv[2]
spanfile = sys.argv[3]
n_decoys = int(sys.argv[4])

os.makedirs(output_dir, exist_ok=True)

pyrosetta.init('-mute all -mp:lipids:has_pore false')

# Fixed histidine pose numbers (haem-coordinating)
fixed_his = [9, 37, 67, 95, 121, 149, 179, 207]

pose = pyrosetta.pose_from_pdb(input_pdb)

from pyrosetta.rosetta.protocols.membrane import AddMembraneMover
amm = AddMembraneMover(spanfile)
amm.apply(pose)

from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
sfxn = ScoreFunctionFactory.create_score_function('franklin2019')

from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    RestrictToRepacking, OperateOnResidueSubset,
    PreventRepackingRLT, RestrictToRepackingRLT
)
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.protocols.relax import FastRelax

# Score storage
scores = []

for i in range(n_decoys):
    work_pose = pose.clone()

    # Set up task factory — fix haem-coordinating HIS, design everything else
    tf = TaskFactory()
    from pyrosetta.rosetta.core.pack.task.operation import InitializeFromCommandline
    tf.push_back(InitializeFromCommandline())

    # Prevent repacking of fixed HIS residues
    fixed_selector = ResidueIndexSelector(','.join(str(r) for r in fixed_his))
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), fixed_selector))

    fr = FastRelax()
    fr.set_scorefxn(sfxn)
    fr.set_task_factory(tf)
    fr.apply(work_pose)

    score = sfxn(work_pose)
    out_path = os.path.join(output_dir, f'decoy_{i:04d}.pdb')
    work_pose.dump_pdb(out_path)
    scores.append((score, out_path))
    print(f'Decoy {i:04d}: score {score:.3f}')
    sys.stdout.flush()

# Write scores file
scores.sort()
with open(os.path.join(output_dir, 'scores.txt'), 'w') as f:
    f.write('rank\tscore\tfile\n')
    for rank, (score, path) in enumerate(scores):
        f.write(f'{rank}\t{score:.3f}\t{os.path.basename(path)}\n')

print(f'Done. Best score: {scores[0][0]:.3f} — {scores[0][1]}')
