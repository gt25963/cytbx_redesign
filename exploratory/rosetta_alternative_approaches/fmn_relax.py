#!/usr/bin/env python3
"""
fmn_relax.py (RQ2)
FastRelax with FMN present, using the validated FMN.params. Backbone is held by coordinate constraints so relaxation opens the pocket locally without deforming the 4-helix bundle. Reports Ca-RMSD to the input.
Usage:
  conda activate pyrosetta
  python fmn_relax.py \
      --in design/FMN_pocket/cycle_1/replace/fmn_final_prerelax.pdb \
      --params ligand_params/FMN.params \
      --out design/FMN_pocket/cycle_1/replace/fmn_final_relaxed.pdb \
      --constrain-bb 0.5 --cycles 3
"""
import argparse

def main():
    #Parse: input structure, ligand params, output, path, backbone constraint strength, relax cycles, and a repack-only mode
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--constrain-bb", type=float, default=0.5)
    ap.add_argument("--cycles", type=int, default=3)
    ap.add_argument("--repack-only", action="store_true")
    args = ap.parse_args()

    import pyrosetta
    from pyrosetta import pose_from_file, get_fa_scorefxn
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.scoring import CA_rmsd, ScoreType

    #Initialise PyRosetta with the FMN ligand params loaded, so FMN is recognised as a valid residue rather than triggering an error
    pyrosetta.init(
        f"-extra_res_fa {args.params} "
        f"-mute all "
        f"-ignore_unrecognized_res false "
        f"-load_PDB_components true"
    )

    #Load the structure and keep an unmodified clone for RMSD comparison later
    pose = pose_from_file(args.inp)
    ref = pose.clone()

    #Standard full-atom score function, with coordinate constraints weighted in
    sfxn = get_fa_scorefxn()
    sfxn.set_weight(ScoreType.coordinate_constraint, 1.0)

    #Add backbone coordinate constraints so relax can open the pocket locally without deforming the overall 4-helix bundle fold
    from pyrosetta.rosetta.protocols.constraint_generator import (
        CoordinateConstraintGenerator, AddConstraints)
    cg = CoordinateConstraintGenerator()
    cg.set_sd(args.constrain_bb)
    cg.set_ca_only(False)
    cg.set_sidechain(False)
    addcst = AddConstraints()
    addcst.add_generator(cg)
    addcst.apply(pose)

    #Two modes: a lighter repack-only pass (side chains only, no backbone movement), or a full constrained FastRelax
    if args.repack_only:
        from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
        from pyrosetta.rosetta.core.pack.task import TaskFactory
        from pyrosetta.rosetta.core.pack.task.operation import (
            RestrictToRepacking, InitializeFromCommandline)
        tf = TaskFactory()
        tf.push_back(InitializeFromCommandline())
        tf.push_back(RestrictToRepacking())
        prm = PackRotamersMover(sfxn)
        prm.task_factory(tf)
        prm.apply(pose)
    else:
        fr = FastRelax(sfxn, args.cycles)
        fr.apply(pose)

    #Save the relaxed structure and report Ca-RMSD plus score change, to confirm the fold was preserved and check the energy improvement
    pose.dump_pdb(args.out)
    rmsd = CA_rmsd(ref, pose)
    e_before = sfxn(ref)
    e_after = sfxn(pose)
    print(f"input:   {args.inp}")
    print(f"output:  {args.out}")
    print(f"Ca-RMSD to input:   {rmsd:.3f} A   (low = fold preserved)")
    print(f"score before/after: {e_before:.1f} / {e_after:.1f} REU")

if __name__ == "__main__":
    main()
