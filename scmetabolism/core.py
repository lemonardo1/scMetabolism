"""
Core functionality for single-cell metabolism analysis.
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler
import warnings
from typing import Optional, Union, Dict, List
import pkg_resources
import gseapy as gp
from concurrent.futures import ThreadPoolExecutor
import os

from .utils import load_gene_sets, alra_imputation


class ScMetabolism:
    """
    Main class for single-cell metabolism analysis.
    
    This class provides methods to quantify metabolic activity at single-cell resolution
    using various enrichment methods including AUCell, ssGSEA, and GSVA-like approaches.
    """
    
    def __init__(self):
        self.metabolism_scores = None
        self.gene_sets = None
        self.method_used = None
        
    def compute_metabolism(
        self,
        count_matrix: Union[pd.DataFrame, np.ndarray],
        method: str = "aucell",
        imputation: bool = False,
        n_cores: int = 2,
        metabolism_type: str = "KEGG"
    ) -> pd.DataFrame:
        """
        Compute metabolism scores for single cells.
        
        Parameters:
        -----------
        count_matrix : pd.DataFrame or np.ndarray
            Gene expression count matrix (genes x cells)
        method : str, default "aucell"
            Method for scoring: "aucell", "ssgsea", "gsva"
        imputation : bool, default False
            Whether to perform ALRA imputation
        n_cores : int, default 2
            Number of cores for parallel processing
        metabolism_type : str, default "KEGG"
            Type of metabolism gene sets: "KEGG" or "REACTOME"
            
        Returns:
        --------
        pd.DataFrame
            Metabolism scores (pathways x cells)
        """
        
        print(f"Your choice is: {metabolism_type}")
        
        # Convert to DataFrame if numpy array
        if isinstance(count_matrix, np.ndarray):
            count_matrix = pd.DataFrame(count_matrix)
            
        # Load gene sets
        self.gene_sets = load_gene_sets(metabolism_type)
        
        # Imputation if requested
        if imputation:
            print("Start imputation...")
            print("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv.")
            count_matrix = alra_imputation(count_matrix)
            
        print("Start quantify the metabolism activity...")
        
        # Compute scores based on method
        if method.lower() == "aucell":
            scores = self._compute_aucell(count_matrix, n_cores)
        elif method.lower() == "ssgsea":
            scores = self._compute_ssgsea(count_matrix, n_cores)
        elif method.lower() == "gsva":
            scores = self._compute_gsva(count_matrix, n_cores)
        else:
            raise ValueError(f"Unknown method: {method}. Choose from 'aucell', 'ssgsea', 'gsva'")
            
        self.metabolism_scores = scores
        self.method_used = method
        
        print("\nPlease Cite:")
        print("Yingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021.")
        print("https://pubmed.ncbi.nlm.nih.gov/34417225/\n")
        
        return scores
    
    def _compute_aucell(self, count_matrix: pd.DataFrame, n_cores: int) -> pd.DataFrame:
        """Compute AUCell scores."""
        from .aucell import AUCell
        
        aucell = AUCell()
        rankings = aucell.build_rankings(count_matrix, n_cores=n_cores)
        scores = aucell.calc_auc(self.gene_sets, rankings)
        
        return scores
    
    def _compute_ssgsea(self, count_matrix: pd.DataFrame, n_cores: int) -> pd.DataFrame:
        """Compute ssGSEA scores."""
        scores_list = []
        
        def compute_pathway_score(pathway_genes, expression_data):
            """Compute ssGSEA score for a single pathway."""
            # Filter genes present in the data
            pathway_genes = [g for g in pathway_genes if g in expression_data.index]
            if len(pathway_genes) < 2:
                return np.zeros(expression_data.shape[1])
                
            # Rank genes for each cell
            ranks = expression_data.rank(axis=0, method='average')
            
            # Calculate enrichment score
            scores = []
            for cell_idx in range(expression_data.shape[1]):
                cell_ranks = ranks.iloc[:, cell_idx].sort_values(ascending=False)
                pathway_ranks = [cell_ranks[g] for g in pathway_genes if g in cell_ranks.index]
                
                if len(pathway_ranks) == 0:
                    scores.append(0)
                    continue
                    
                # Simplified ssGSEA calculation
                n_genes = len(cell_ranks)
                pathway_positions = sorted([n_genes - r + 1 for r in pathway_ranks])
                
                # Calculate enrichment score
                es = 0
                for i, pos in enumerate(pathway_positions):
                    es += (pos - i * n_genes / len(pathway_positions)) / n_genes
                    
                scores.append(es)
                
            return np.array(scores)
        
        # Process pathways
        with ThreadPoolExecutor(max_workers=n_cores) as executor:
            futures = {}
            for pathway_name, pathway_genes in self.gene_sets.items():
                future = executor.submit(compute_pathway_score, pathway_genes, count_matrix)
                futures[pathway_name] = future
                
            for pathway_name, future in futures.items():
                scores_list.append((pathway_name, future.result()))
        
        # Create DataFrame
        scores_dict = {name: scores for name, scores in scores_list}
        scores_df = pd.DataFrame(scores_dict, index=count_matrix.columns).T
        
        return scores_df
    
    def _compute_gsva(self, count_matrix: pd.DataFrame, n_cores: int) -> pd.DataFrame:
        """Compute GSVA-like scores."""
        # Simplified GSVA implementation
        # Normalize data
        normalized_data = count_matrix.div(count_matrix.sum(axis=0), axis=1) * 1e6
        normalized_data = np.log2(normalized_data + 1)
        
        scores_dict = {}
        
        for pathway_name, pathway_genes in self.gene_sets.items():
            # Filter genes present in data
            pathway_genes = [g for g in pathway_genes if g in normalized_data.index]
            if len(pathway_genes) < 2:
                scores_dict[pathway_name] = np.zeros(normalized_data.shape[1])
                continue
                
            # Get pathway expression
            pathway_expr = normalized_data.loc[pathway_genes]
            
            # Calculate scores (simplified GSVA)
            pathway_scores = []
            for cell_idx in range(normalized_data.shape[1]):
                cell_expr = normalized_data.iloc[:, cell_idx]
                pathway_cell_expr = pathway_expr.iloc[:, cell_idx]
                
                # Rank-based enrichment
                all_ranks = cell_expr.rank(pct=True)
                pathway_ranks = all_ranks[pathway_genes]
                
                # Calculate enrichment score
                score = pathway_ranks.mean() - 0.5
                pathway_scores.append(score)
                
            scores_dict[pathway_name] = np.array(pathway_scores)
        
        scores_df = pd.DataFrame(scores_dict, index=count_matrix.columns).T
        return scores_df
    
    def compute_metabolism_scanpy(
        self,
        adata,
        method: str = "aucell",
        imputation: bool = False,
        n_cores: int = 2,
        metabolism_type: str = "KEGG",
        layer: Optional[str] = None
    ):
        """
        Compute metabolism scores for AnnData object (scanpy compatible).
        
        Parameters:
        -----------
        adata : AnnData
            Annotated data object
        method : str, default "aucell"
            Method for scoring
        imputation : bool, default False
            Whether to perform imputation
        n_cores : int, default 2
            Number of cores
        metabolism_type : str, default "KEGG"
            Gene set type
        layer : str, optional
            Layer to use for computation
            
        Returns:
        --------
        AnnData
            Updated AnnData object with metabolism scores in .obsm['metabolism']
        """
        
        # Extract count matrix
        if layer is None:
            count_matrix = pd.DataFrame(
                adata.X.T if hasattr(adata.X, 'todense') else adata.X.T,
                index=adata.var_names,
                columns=adata.obs_names
            )
        else:
            count_matrix = pd.DataFrame(
                adata.layers[layer].T,
                index=adata.var_names,
                columns=adata.obs_names
            )
        
        # Compute scores
        scores = self.compute_metabolism(
            count_matrix, method, imputation, n_cores, metabolism_type
        )
        
        # Add to AnnData
        adata.obsm['metabolism'] = scores.T
        adata.uns['metabolism_pathways'] = list(scores.index)
        
        return adata