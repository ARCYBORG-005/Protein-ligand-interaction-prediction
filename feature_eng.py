from rdkit.Chem import AllChem
import random

class MoleculeFeaturizer:
    @staticmethod
    def smiles_augmentation(smiles):
        """Randomizes SMILES string to improve GNN robustness."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return smiles
        return Chem.MolToSmiles(mol, doRandom=True, canonical=False)

    @staticmethod
    def get_morgan_fingerprint(smiles, radius=2, n_bits=2048):
        """Generates Morgan Fingerprints as bit vectors."""
        mol = Chem.MolFromSmiles(smiles)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        return np.array(fp)
