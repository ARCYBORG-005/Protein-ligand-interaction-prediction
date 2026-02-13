import pandas as pd
from rdkit import Chem
from scipy import stats
import numpy as np

class DataCleaner:
    def __init__(self, affinity_col='pkd'):
        self.affinity_col = affinity_col

    def remove_invalid_smiles(self, df):
        """Removes rows with chemically invalid SMILES strings."""
        return df[df['smiles'].apply(lambda x: Chem.MolFromSmiles(x) is not None)]

    def handle_outliers(self, df, threshold=3):
        """Uses Z-score to remove extreme binding affinity outliers."""
        z_scores = np.abs(stats.zscore(df[self.affinity_col]))
        return df[z_scores < threshold]

    def standardize_protein(self, df):
        """Removes non-standard amino acids and short sequences."""
        valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
        df['protein_sequence'] = df['protein_sequence'].str.upper()
        mask = df['protein_sequence'].apply(lambda x: all(c in valid_aas for c in x) and len(x) > 50)
        return df[mask]
