"""
TRAINING AND EVALUATION OUTPUT (PDBbind Dataset v2020)
------------------------------------------------------
[Phase 1: Data Cleaning]
- Initial SMILES count: 19,443
- Invalid SMILES removed: 142
- Binding affinity outliers (Z > 3.0) removed: 86
- Final training set size: 19,215

[Phase 2: Training Progress]
Epoch [01/100] | Loss: 2.145 | MSE: 2.145 | R2: -0.12 | Time: 45s
Epoch [25/100] | Loss: 0.892 | MSE: 0.890 | R2: 0.58  | Time: 42s
Epoch [50/100] | Loss: 0.541 | MSE: 0.538 | R2: 0.71  | Time: 43s (Learning Rate Reduced)
Epoch [100/100]| Loss: 0.421 | MSE: 0.419 | R2: 0.79  | Time: 42s

[Phase 3: Final Test Set Metrics]
- Mean Squared Error (MSE): 0.4215 (Traditional Docking: 1.84)
- R-Squared (R2):           0.7854 (Traditional Docking: 0.45)
- Pearson Correlation:      0.8922
- Concordance Index (CI):   0.8741

[Phase 4: Interpretability Output]
- GNNExplainer: Top 3 atomic contributors identified for Ligand ID #402:
  1. Carbonyl Oxygen (Node 4) - Weight: 0.88
  2. Phenyl Ring C2 (Node 12) - Weight: 0.74
  3. Hydroxyl Group (Node 8) - Weight: 0.69
"""
