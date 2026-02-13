import torch.nn as nn
from torch_geometric.nn import GCNConv, global_mean_pool

class HybridDrugTargetNet(nn.Module):
    def __init__(self):
        super().__init__()
        # GNN for Ligand
        self.ligand_gnn = nn.ModuleList([
            GCNConv(78, 128),
            GCNConv(128, 256)
        ])
        
        # CNN for Protein
        self.protein_cnn = nn.Sequential(
            nn.Embedding(25, 128),
            nn.Conv1d(128, 64, kernel_size=8),
            nn.ReLU(),
            nn.AdaptiveMaxPool1d(1)
        )

        # Cross-Attention / Fusion
        self.fc_layers = nn.Sequential(
            nn.Linear(256 + 64, 512),
            nn.Dropout(0.2),
            nn.ReLU(),
            nn.Linear(512, 1) # Binding Affinity pKd
        )

    def forward(self, ligand_data, protein_seq):
        # 1. Ligand Branch
        x, edge_index, batch = ligand_data.x, ligand_data.edge_index, ligand_data.batch
        for conv in self.ligand_gnn:
            x = torch.relu(conv(x, edge_index))
        ligand_feat = global_mean_pool(x, batch)

        # 2. Protein Branch
        protein_feat = self.protein_cnn(protein_seq).squeeze(-1)

        # 3. Concatenate and Predict
        combined = torch.cat([ligand_feat, protein_feat], dim=1)
        return self.fc_layers(combined)
