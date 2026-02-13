from torch_geometric.explain import Explainer, GNNExplainer

def explain_prediction(model, drug_sample):
    """Explains which atoms in the ligand are most important for binding."""
    explainer = Explainer(
        model=model,
        algorithm=GNNExplainer(epochs=200),
        explanation_type='model',
        node_mask_type='attributes',
        edge_mask_type='object',
        model_config=dict(mode='regression', task_level='graph', return_type='raw'),
    )
    explanation = explainer(drug_sample.x, drug_sample.edge_index)
    return explanation.node_mask # Importance scores per atom
