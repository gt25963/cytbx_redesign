"""
Generate Rosetta .params files for FMN and U10 (ubiquinone-10) for RQ2.
SMILES sourced and verified from RCSB PDB Chemical Component Dictionary.
"""

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
    if embed_result != 0:
        print(f"WARNING: initial embedding failed for {name}, retrying with random coords")
        embed_result = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
    AllChem.MMFFOptimizeMolecule(mol)

    p = Params.from_mol(mol, name=name)
    p.dump(f"{name}.params")
    p.dump_pdb_conf(f"{name}_conf.pdb")
    print(f"Wrote {name}.params and {name}_conf.pdb")

print("Done.")
