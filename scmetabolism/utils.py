"""
Utility functions for scMetabolism package.
"""

import numpy as np
import pandas as pd
from scipy.sparse.linalg import svds
from sklearn.decomposition import TruncatedSVD
import pkg_resources
import os
from typing import Dict, List, Union


def load_gene_sets(metabolism_type: str = "KEGG") -> Dict[str, List[str]]:
    """
    Load gene sets from various sources.
    
    Parameters:
    -----------
    metabolism_type : str
        Type of gene sets: 
        - "KEGG": KEGG metabolism pathways
        - "REACTOME": REACTOME metabolism pathways  
        - "GO_metabolism": GO metabolism-related terms
        - "GO_all": All GO terms (BP, MF, CC)
        - "GO_BP": GO Biological Process
        - "GO_MF": GO Molecular Function
        - "GO_CC": GO Cellular Component
        
    Returns:
    --------
    Dict[str, List[str]]
        Dictionary mapping pathway names to gene lists
    """
    
    # Handle traditional GMT files (KEGG, REACTOME)
    if metabolism_type.upper() in ["KEGG", "REACTOME"]:
        if metabolism_type.upper() == "KEGG":
            gmt_file = "KEGG_metabolism_nc.gmt"
        else:
            gmt_file = "REACTOME_metabolism.gmt"
        
        # Try to load from package data
        try:
            gmt_path = pkg_resources.resource_filename('scmetabolism', f'data/{gmt_file}')
        except:
            # Fallback to local data directory
            gmt_path = os.path.join('data', gmt_file)
        
        if not os.path.exists(gmt_path):
            raise FileNotFoundError(f"Gene set file not found: {gmt_path}")
        
        gene_sets = {}
        
        with open(gmt_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    pathway_name = parts[0]
                    genes = parts[2:]  # Skip description field
                    gene_sets[pathway_name] = genes
        
        return gene_sets
    
    # Handle GO gene sets
    elif metabolism_type.startswith("GO_"):
        from .go_analysis import download_go_gene_sets
        
        # Map metabolism_type to GO analysis parameters
        go_type_map = {
            "GO_metabolism": "metabolism",
            "GO_all": "all", 
            "GO_BP": "biological_process",
            "GO_MF": "molecular_function",
            "GO_CC": "cellular_component"
        }
        
        if metabolism_type not in go_type_map:
            raise ValueError(f"Unknown GO type: {metabolism_type}")
        
        print(f"Downloading GO gene sets: {metabolism_type}")
        gene_sets = download_go_gene_sets(
            organism="human",
            gene_set_type=go_type_map[metabolism_type]
        )
        
        return gene_sets
    
    else:
        raise ValueError(
            f"Unknown metabolism_type: {metabolism_type}. "
            f"Choose from: KEGG, REACTOME, GO_metabolism, GO_all, GO_BP, GO_MF, GO_CC"
        )


def alra_imputation(count_matrix: pd.DataFrame, k: int = None) -> pd.DataFrame:
    """
    ALRA (Adaptively-thresholded Low Rank Approximation) imputation.
    
    Parameters:
    -----------
    count_matrix : pd.DataFrame
        Gene expression matrix (genes x cells)
    k : int, optional
        Number of components for SVD. If None, automatically determined.
        
    Returns:
    --------
    pd.DataFrame
        Imputed expression matrix
    """
    
    # Convert to numpy array
    X = count_matrix.values.astype(float)
    
    # Log transform (add pseudocount)
    X_log = np.log2(X + 1)
    
    # Center the data
    X_centered = X_log - np.mean(X_log, axis=1, keepdims=True)
    
    # Determine k if not provided
    if k is None:
        # Use elbow method to find optimal k
        k = min(100, min(X.shape) - 1)
        
        # Perform SVD with different k values to find elbow
        explained_var_ratios = []
        k_values = range(10, min(k, 50), 5)
        
        for test_k in k_values:
            try:
                svd = TruncatedSVD(n_components=test_k, random_state=42)
                svd.fit(X_centered.T)
                explained_var_ratios.append(np.sum(svd.explained_variance_ratio_))
            except:
                break
        
        # Find elbow point
        if len(explained_var_ratios) > 2:
            diffs = np.diff(explained_var_ratios)
            k = k_values[np.argmin(diffs[1:]) + 1]
        else:
            k = 20  # Default fallback
    
    # Perform SVD
    try:
        U, s, Vt = svds(X_centered.T, k=k)
        
        # Reconstruct low-rank approximation
        X_lr = (U * s) @ Vt
        X_lr = X_lr.T
        
        # Add back the mean
        X_lr = X_lr + np.mean(X_log, axis=1, keepdims=True)
        
        # Convert back from log space
        X_imputed = np.power(2, X_lr) - 1
        
        # Preserve zeros where original data was zero
        mask = (X == 0)
        X_imputed[mask] = 0
        
        # Create DataFrame with original index and columns
        result = pd.DataFrame(
            X_imputed,
            index=count_matrix.index,
            columns=count_matrix.columns
        )
        
        return result
        
    except Exception as e:
        print(f"ALRA imputation failed: {e}")
        print("Returning original matrix")
        return count_matrix


def preprocess_data(
    count_matrix: pd.DataFrame,
    min_genes: int = 200,
    min_cells: int = 3,
    normalize: bool = True,
    log_transform: bool = True
) -> pd.DataFrame:
    """
    Basic preprocessing of single-cell data.
    
    Parameters:
    -----------
    count_matrix : pd.DataFrame
        Raw count matrix (genes x cells)
    min_genes : int
        Minimum number of genes per cell
    min_cells : int
        Minimum number of cells per gene
    normalize : bool
        Whether to normalize to counts per million
    log_transform : bool
        Whether to log transform
        
    Returns:
    --------
    pd.DataFrame
        Preprocessed matrix
    """
    
    # Filter cells with too few genes
    genes_per_cell = (count_matrix > 0).sum(axis=0)
    valid_cells = genes_per_cell >= min_genes
    count_matrix = count_matrix.loc[:, valid_cells]
    
    # Filter genes expressed in too few cells
    cells_per_gene = (count_matrix > 0).sum(axis=1)
    valid_genes = cells_per_gene >= min_cells
    count_matrix = count_matrix.loc[valid_genes, :]
    
    if normalize:
        # Normalize to counts per million
        count_matrix = count_matrix.div(count_matrix.sum(axis=0), axis=1) * 1e6
    
    if log_transform:
        # Log transform with pseudocount
        count_matrix = np.log2(count_matrix + 1)
    
    return count_matrix


def calculate_pathway_activity(
    expression_data: pd.DataFrame,
    gene_sets: Dict[str, List[str]],
    method: str = "mean"
) -> pd.DataFrame:
    """
    Calculate pathway activity scores.
    
    Parameters:
    -----------
    expression_data : pd.DataFrame
        Gene expression data (genes x cells)
    gene_sets : Dict[str, List[str]]
        Gene sets dictionary
    method : str
        Method for aggregation: "mean", "median", "sum"
        
    Returns:
    --------
    pd.DataFrame
        Pathway activity scores (pathways x cells)
    """
    
    pathway_scores = {}
    
    for pathway_name, genes in gene_sets.items():
        # Find genes present in the data
        present_genes = [g for g in genes if g in expression_data.index]
        
        if len(present_genes) == 0:
            # No genes found, set score to 0
            pathway_scores[pathway_name] = np.zeros(expression_data.shape[1])
        else:
            # Calculate pathway score
            pathway_expr = expression_data.loc[present_genes]
            
            if method == "mean":
                scores = pathway_expr.mean(axis=0)
            elif method == "median":
                scores = pathway_expr.median(axis=0)
            elif method == "sum":
                scores = pathway_expr.sum(axis=0)
            else:
                raise ValueError("method must be 'mean', 'median', or 'sum'")
                
            pathway_scores[pathway_name] = scores.values
    
    # Create DataFrame
    scores_df = pd.DataFrame(
        pathway_scores,
        index=expression_data.columns
    ).T
    
    return scores_df