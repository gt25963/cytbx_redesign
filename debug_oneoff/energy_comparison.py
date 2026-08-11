import pyrosetta
import sys
pyrosetta.init('-mute all -mp:lipids:has_pore false')

#Diagnostic comparison: relax two C2 dimer designs that differ at position 75 (Tyr vs Phe) to directly compare per-residue energy at that position, testing whether Tyr75 alone explains the observed score gap between high- and low-scoring C2 sequences (Results 1.2, "Position 75" finding)
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

#Membrane-aware Rosetta score function used throughout RQ1's Rosetta trials
sfxn = ScoreFunctionFactory.create_score_function('franklin2019')

#Two representative structures: one with Tyr75 (high-scoring pattern), one with Phe75 (regressed pattern)
structures = {
    "id16_Y75": "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/diagnostics/energy_comparison/id16_Y75_holo.pdb",
    "id41_F75": "/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline/diagnostics/energy_comparison/id41_F75_holo.pdb",
}

results = {}
for name, path in structures.items():
    print(f"--- Processing {name} ---")
    pose = pyrosetta.pose_from_pdb(path)

    #Set up membrane environment from the spanfile before scoring
    amm = AddMembraneMover(spanfile)
    amm.apply(pose)

    #Hold haem-coordinating histidines fixed during relax, same as the main design pipeline
    tf = TaskFactory()
    fixed_selector = ResidueIndexSelector(','.join(str(r) for r in fixed_his))
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), fixed_selector))

    #Relax the structure under the membrane score function
    fr = FastRelax()
    fr.set_scorefxn(sfxn)
    fr.set_task_factory(tf)
    fr.apply(pose)

    total_score = sfxn(pose)
    print(f"{name} total score: {total_score:.3f}")

    #Pull per-residue energy specifically at position 75 on both chains (C2 = 2 copies)
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

#Print side-by-side comparison of total and position-75 energies for both variants
print("\n=== SUMMARY ===")
for name, r in results.items():
    print(f"{name}: total={r['total_score']:.3f}, res75_A={r['res75_A']:.3f}, res75_B={r['res75_B']:.3f}")
