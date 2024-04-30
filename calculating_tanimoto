import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def main():

    smiles1 = "CCCCOCOC(=O)NC"
    smiles2 = "CNC(=O)OCCCCC(C)C"
    smiles3 = "c1ccccc1OC(=O)O"


    mol = [Chem.MolFromSmiles(smiles1), Chem.MolFromSmiles(smiles2), Chem.MolFromSmiles(smiles3)]
    fp = [AllChem.GetMorganFingerprintAsBitVect(x, 3, nBits=2048) for x in mol]

    print(f"Tanimoto similarity between molecules I  & II  : {calc_fp(fp[0], fp[1])}", flush=True)
    print(f"Tanimoto similarity between molecules I  & III : {calc_fp(fp[0], fp[2])}", flush=True)
    print(f"Tanimoto similarity between molecules II & III : {calc_fp(fp[1], fp[2])}", flush=True)


def calc_fp(fp_a, fp_b):

    # fp_a: fingerprint of molecule a
    # fp_b: fingerprint of molecule b

    ###################################
    ######## Change this part #########
    ###################################

    Tsim = 0.0 
    top = np.dot(fp_a, fp_b)
    bot = sum(fp_a) + sum(fp_b) - top
    Tsim = top / bot
    return Tsim

if __name__ == "__main__":
    main()


