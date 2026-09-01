#Generate Rosetta .params files for FMN and U10 (ubiquinone-10) for RQ2. - work prior to the actual Q8 cofactor
#SMILES sourced and verified from RCSB PDB Chemical Component Dictionary.

# IMPORTANT: rdkit_to_params incorrectly applies ring-conformer treatment (ADD_RING + PROPERTIES CYCLIC) to aromatic/rigid rings
# This affects both FMN's isoalloxazine and U10/Q8's quinone ring.
# The .params files this script writes are NOT yet PyRosetta-loadable as-is; 
## The three offending lines must be stripped (via sed) as a separate post-processing step before use. 
## The validated, corrected FMN.params/Q8.params actually used in the pipeline are in rq2/ligand_params/ - if re-running this script, don't skip that step.

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_to_params import Params

ligands = {
    "FMN": "Cc1cc2c(cc1C)N(C3=NC(=O)NC(=O)C3=N2)C[C@@H]([C@@H]([C@@H](COP(=O)(O)O)O)O)O",
    "U10": "CC1=C(C(=O)C(=C(C1=O)OC)OC)C\\C=C(/C)\\CC\\C=C(/C)\\CC\\C=C(/C)\\CC\\C=C(/C)\\CC\\C=C(/C)\\CC\\C=C(/C)\\CC\\C=C(/C)\\CC\\C=C(/C)\\CCC=C(C)C",
}

for name, smiles in ligands.items():
    print(f"--- Processing {name} ---")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"ERROR: RDKit could not parse SMILES for {name}")
        continue

    mol = AllChem.AddHs(mol)
    embed_result = AllChem.EmbedMolecule(mol, randomSeed=42)
    ## random-coords fallback = standard RDKit workaround
    if embed_result != 0:
        print(f"WARNING: initial embedding failed for {name}, retrying with random coords")
        embed_result = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
    AllChem.MMFFOptimizeMolecule(mol)

    p = Params.from_mol(mol, name=name)
    p.dump(f"{name}.params") ## needs the ring-conformer sed fix before use
    p.dump_pdb_conf(f"{name}_conf.pdb")
    print(f"Wrote {name}.params and {name}_conf.pdb")

print("Done.")
