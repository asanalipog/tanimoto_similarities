import numpy as np

from rdkit import Chem
from rdkit.Chem import Draw, AllChem

def main():

    smiles1 = "CCCCOCOC(=O)NC"
    smiles2 = "CNC(=O)OCCCCC(C)C"
    smiles3 = "c1ccccc1OC(=O)O"

    mol = Chem.MolFromSmiles(smiles1)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=2048)
    print(f"Fingerprint of molecule #1")
    print([fp[x] for x in range(2048)])
    fp_image(smiles1, "smiles1")

    mol = Chem.MolFromSmiles(smiles2)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=2048)
    print(f"\nFingerprint of molecule #2")
    print([fp[x] for x in range(2048)])
    fp_image(smiles2, "smiles2")

    mol = Chem.MolFromSmiles(smiles3)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=2048)
    print(f"\nFingerprint of molecule #3")
    print([fp[x] for x in range(2048)])
    fp_image(smiles3, "smiles3")


def fp_image(smiles, filename):
    mol = Chem.MolFromSmiles(smiles)
    bi = {}
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, bitInfo=bi)
    tpls = [(mol,x,bi) for x in fp.GetOnBits()]
    mfp2_svg = Draw.DrawMorganBits(tpls[:],molsPerRow=6,legends=[str(x) for x in fp.GetOnBits()][:])
    with open(f"{filename}.svg", "w") as f:
        f.write(mfp2_svg)

if __name__ == "__main__":
    main()

