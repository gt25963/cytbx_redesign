import pyrosetta
import sys

pyrosetta.init('-mute all -mp:lipids:has_pore false')

spanfile = "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/cytbx_C2_top2_dimer.span"
fixed_his = [9, 37, 67, 95, 121, 149, 179, 207]

from pyrosetta.rosetta.protocols.membrane import AddMembraneMover
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    OperateOnResidueSubset, PreventRepackingRLT
)
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.protocols.relax import FastRelax

sfxn = ScoreFunctionFactory.create_score_function('franklin2019')

structures = {
    "id16_Y75": "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/diagnostics/energy_comparison/id16_Y75_holo.pdb",
    "id41_F75": "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/diagnostics/energy_comparison/id41_F75_holo.pdb",
}

results = {}

for name, path in structures.items():
    print(f"--- Processing {name} ---")
    pose = pyrosetta.pose_from_pdb(path)
    amm = AddMembraneMover(spanfile)
    amm.apply(pose)

    tf = TaskFactory()
    fixed_selector = ResidueIndexSelector(','.join(str(r) for r in fixed_his))
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), fixed_selector))

    fr = FastRelax()
    fr.set_scorefxn(sfxn)
    fr.set_task_factory(tf)
    fr.apply(pose)

    total_score = sfxn(pose)
    print(f"{name} total score: {total_score:.3f}")

    energies = pose.energies()
    res75_A_energy = energies.residue_total_energy(75)
    chain_a_length = pose.chain_end(1)
    res75_B_pose_num = chain_a_length + 75
    res75_B_energy = energies.residue_total_energy(res75_B_pose_num)

    print(f"{name} - Position 75 (chain A) per-residue energy: {res75_A_energy:.3f}")
    print(f"{name} - Position 75 (chain B) per-residue energy: {res75_B_energy:.3f}")

    pose.dump_pdb(f"/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/diagnostics/energy_comparison/{name}_relaxed.pdb")

    results[name] = {
        "total_score": total_score,
        "res75_A": res75_A_energy,
        "res75_B": res75_B_energy,
    }

print("\n=== SUMMARY ===")
for name, r in results.items():
    print(f"{name}: total={r['total_score']:.3f}, res75_A={r['res75_A']:.3f}, res75_B={r['res75_B']:.3f}")
