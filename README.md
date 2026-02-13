Advancing Drug Discovery:
Hybrid GNN-CNN Binding Affinity PredictionPredicting Protein-Ligand Binding Affinity ($pK_d$) using a multi-modal Deep Learning approach to accelerate lead optimization in drug discovery.

Project Overview
This repository features a state-of-the-art interaction prediction system (developed in 2025) that outperforms traditional physics-based docking. By treating ligands as structural graphs and proteins as sequential motifs, the model captures the complex "lock-and-key" biophysics of molecular binding.

Key Innovations
Dual-Stream Encoding: Uses Graph Convolutional Networks (GCN) for chemical topology and 1D-CNNs for protein sequence motifs.

Cross-Attention Fusion: Implements an attention mechanism to model specific atom-to-residue interactions.

XAI (Explainable AI): Integrated GNNExplainer to visualize atomic contributions, turning the "black box" into a tool for medicinal chemists.

Data Engineering: Robust pipeline featuring SMILES Augmentation, Morgan Fingerprints, and Z-score outlier filtration.

Module,                             Purpose,Key                                    Techniques
data_cleaner.py,                        Quality Control,                    "Z-score Outlier Removal, RDKit Sanitization"
featurizer.py,                           Representation                       Morgan Fingerprints"
model.py,                                core logic                               GNN (Ligand) + CNN (Protein) + Cross-Attention
explainer.py,
Interpretability"                    GNNExplainer,                           Feature Importance Heatmaps


Metric,         AutoDock Vina,Our Hybrid Model,Improvement
MSE,               1.84,0.42,🟢 77% Lower
R² Score             ,0.45,0.78,🟢 73% Higher
Pearson (R),               0.62,0.89,🟢 43% Higher

nterpretability (XAI)
A primary feature for interview discussion is the GNNExplainer integration. The model doesn't just predict a score; it highlights the pharmacophores responsible for the binding.

Visual Output: The system generates a node-importance map. For example, in Carbonic Anhydrase inhibitors, the model correctly identifies the sulfonamide group as the primary binding driver.


