from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdMolDescriptors

fmn_smiles = "CC1=CC2=C(C=C1C)N(C3=NC(=O)NC(=O)C3=N2)CC(C(C(COP(=O)(O)O)O)O)O"
q8_smiles = "CC1=C(C(=O)C(=C(C1=O)OC)OC)CC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C"

def render(smiles, name, width, height):
    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(width, height)
    d.drawOptions().bondLineWidth = 2
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    with open(f"{name}.png", "wb") as f:
        f.write(d.GetDrawingText())
    print(name, rdMolDescriptors.CalcMolFormula(mol))

render(fmn_smiles, "FMN", 850, 620)
render(q8_smiles, "Q8", 1350, 320)
