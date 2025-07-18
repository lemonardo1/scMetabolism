"""
Basic usage example for scMetabolism package.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scmetabolism import ScMetabolism, MetabolismPlotter

def create_example_data():
    """Create example single-cell RNA-seq data."""
    
    # Create synthetic data
    n_genes = 2000
    n_cells = 500
    
    # Gene names (some real metabolism-related genes)
    metabolism_genes = [
        'ALDOA', 'ENO1', 'GAPDH', 'GPI', 'HK1', 'HK2', 'LDHA', 'PFKL', 'PFKM', 'PFKP',
        'PGK1', 'PGAM1', 'PKM', 'TPI1', 'ACLY', 'FASN', 'ACC1', 'SCD', 'ELOVL6',
        'COX1', 'COX2', 'COX3', 'NDUFA1', 'NDUFB1', 'SDHB', 'UQCRC1', 'ATP5A1'
    ]
    
    other_genes = [f'GENE_{i}' for i in range(len(metabolism_genes), n_genes)]
    all_genes = metabolism_genes + other_genes
    
    # Cell names
    cell_names = [f'Cell_{i}' for i in range(n_cells)]
    
    # Create count matrix with some structure
    np.random.seed(42)
    
    # Base expression levels
    base_expression = np.random.negative_binomial(5, 0.3, size=(n_genes, n_cells))
    
    # Add some cell type structure
    cell_types = ['Type_A', 'Type_B', 'Type_C']
    cells_per_type = n_cells // len(cell_types)
    
    # Boost metabolism genes in different cell types
    for i, cell_type in enumerate(cell_types):
        start_idx = i * cells_per_type
        end_idx = (i + 1) * cells_per_type if i < len(cell_types) - 1 else n_cells
        
        # Boost different metabolism pathways in different cell types
        if cell_type == 'Type_A':  # Glycolysis
            glycolysis_indices = [j for j, gene in enumerate(all_genes) 
                                if gene in ['ALDOA', 'ENO1', 'GAPDH', 'HK1', 'LDHA', 'PKM']]
            base_expression[glycolysis_indices, start_idx:end_idx] *= 3
            
        elif cell_type == 'Type_B':  # Oxidative phosphorylation
            oxphos_indices = [j for j, gene in enumerate(all_genes)
                            if gene in ['COX1', 'COX2', 'NDUFA1', 'SDHB', 'ATP5A1']]
            base_expression[oxphos_indices, start_idx:end_idx] *= 4
            
        elif cell_type == 'Type_C':  # Fatty acid synthesis
            fas_indices = [j for j, gene in enumerate(all_genes)
                         if gene in ['ACLY', 'FASN', 'ACC1', 'SCD']]
            base_expression[fas_indices, start_idx:end_idx] *= 2.5
    
    # Create DataFrame
    count_matrix = pd.DataFrame(
        base_expression,
        index=all_genes,
        columns=cell_names
    )
    
    # Create metadata
    metadata = pd.DataFrame({
        'cell_type': [cell_types[i // cells_per_type] for i in range(n_cells)],
        'batch': np.random.choice(['Batch1', 'Batch2'], n_cells),
        'n_genes': (count_matrix > 0).sum(axis=0).values
    }, index=cell_names)
    
    return count_matrix, metadata


def main():
    """Main example workflow."""
    
    print("=== scMetabolism Python Package Example ===\n")
    
    # 1. Create example data
    print("1. Creating example data...")
    count_matrix, metadata = create_example_data()
    print(f"   Count matrix shape: {count_matrix.shape}")
    print(f"   Cell types: {metadata['cell_type'].value_counts().to_dict()}")
    
    # 2. Initialize ScMetabolism
    print("\n2. Initializing ScMetabolism...")
    sc_metab = ScMetabolism()
    
    # 3. Compute metabolism scores
    print("\n3. Computing metabolism scores...")
    metabolism_scores = sc_metab.compute_metabolism(
        count_matrix=count_matrix,
        method="aucell",
        imputation=False,
        n_cores=2,
        metabolism_type="KEGG"
    )
    
    print(f"   Computed scores for {metabolism_scores.shape[0]} pathways")
    print(f"   Top pathways by mean activity:")
    mean_activity = metabolism_scores.mean(axis=1).sort_values(ascending=False)
    for pathway in mean_activity.head(5).index:
        print(f"     - {pathway}: {mean_activity[pathway]:.3f}")
    
    # 4. Create visualizations
    print("\n4. Creating visualizations...")
    plotter = MetabolismPlotter()
    
    # Create synthetic 2D embedding for visualization
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    
    # PCA for initial dimensionality reduction
    pca = PCA(n_components=50, random_state=42)
    pca_coords = pca.fit_transform(count_matrix.T)
    
    # t-SNE for visualization
    tsne = TSNE(n_components=2, random_state=42, perplexity=30)
    tsne_coords = tsne.fit_transform(pca_coords)
    
    embedding = pd.DataFrame(
        tsne_coords,
        columns=['tSNE1', 'tSNE2'],
        index=count_matrix.columns
    )
    
    # Box plot
    print("   Creating box plot...")
    pathways_to_plot = mean_activity.head(6).index.tolist()
    
    fig = plotter.box_plot(
        metabolism_scores=metabolism_scores,
        metadata=metadata,
        pathways=pathways_to_plot,
        group_by='cell_type',
        ncols=2,
        figsize=(12, 8)
    )
    plt.savefig('metabolism_boxplot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Dot plot
    print("   Creating dot plot...")
    fig = plotter.dot_plot(
        metabolism_scores=metabolism_scores,
        metadata=metadata,
        pathways=pathways_to_plot[:4],
        group_by='cell_type',
        normalize='row',
        figsize=(8, 6)
    )
    plt.savefig('metabolism_dotplot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Dimension plot
    print("   Creating dimension plot...")
    fig = plotter.dim_plot(
        embedding=embedding,
        metabolism_scores=metabolism_scores,
        pathway=pathways_to_plot[0],
        embedding_type='tsne',
        size=1.5,
        figsize=(8, 6)
    )
    plt.savefig('metabolism_dimplot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Heatmap
    print("   Creating heatmap...")
    fig = plotter.heatmap(
        metabolism_scores=metabolism_scores,
        pathways=pathways_to_plot,
        cluster_rows=True,
        cluster_cols=True,
        figsize=(10, 6)
    )
    plt.savefig('metabolism_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 5. Summary statistics
    print("\n5. Summary statistics:")
    print(f"   Total pathways analyzed: {metabolism_scores.shape[0]}")
    print(f"   Total cells analyzed: {metabolism_scores.shape[1]}")
    print(f"   Mean score range: {metabolism_scores.values.min():.3f} - {metabolism_scores.values.max():.3f}")
    
    # Pathway activity by cell type
    print("\n   Pathway activity by cell type:")
    for cell_type in metadata['cell_type'].unique():
        cell_type_cells = metadata[metadata['cell_type'] == cell_type].index
        cell_type_scores = metabolism_scores.loc[:, cell_type_cells].mean(axis=1)
        top_pathway = cell_type_scores.idxmax()
        print(f"     {cell_type}: Top pathway = {top_pathway} (score: {cell_type_scores[top_pathway]:.3f})")
    
    print("\n=== Analysis complete! ===")
    print("Generated files:")
    print("  - metabolism_boxplot.png")
    print("  - metabolism_dotplot.png") 
    print("  - metabolism_dimplot.png")
    print("  - metabolism_heatmap.png")


if __name__ == "__main__":
    main()