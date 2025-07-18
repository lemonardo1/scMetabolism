"""
AUCell implementation for Python.
"""

import numpy as np
import pandas as pd
from scipy import stats
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Union


class AUCell:
    """
    Python implementation of AUCell algorithm for gene set enrichment analysis.
    
    Based on the R package AUCell by Sara Aibar et al.
    """
    
    def __init__(self):
        self.rankings = None
        
    def build_rankings(
        self,
        expression_matrix: pd.DataFrame,
        n_cores: int = 1,
        plot_stats: bool = False
    ) -> pd.DataFrame:
        """
        Build gene rankings for each cell.
        
        Parameters:
        -----------
        expression_matrix : pd.DataFrame
            Gene expression matrix (genes x cells)
        n_cores : int
            Number of cores for parallel processing
        plot_stats : bool
            Whether to plot ranking statistics
            
        Returns:
        --------
        pd.DataFrame
            Gene rankings for each cell
        """
        
        print("Building gene rankings...")
        
        # Rank genes for each cell (higher expression = higher rank)
        rankings = expression_matrix.rank(axis=0, method='average', ascending=False)
        
        self.rankings = rankings
        
        if plot_stats:
            self._plot_ranking_stats(rankings)
            
        return rankings
    
    def calc_auc(
        self,
        gene_sets: Dict[str, List[str]],
        rankings: pd.DataFrame = None,
        auc_max_rank: Union[int, float] = 0.05
    ) -> pd.DataFrame:
        """
        Calculate AUC scores for gene sets.
        
        Parameters:
        -----------
        gene_sets : Dict[str, List[str]]
            Dictionary of gene sets
        rankings : pd.DataFrame, optional
            Gene rankings. If None, uses stored rankings.
        auc_max_rank : int or float
            Maximum rank to consider for AUC calculation.
            If float, interpreted as fraction of total genes.
            
        Returns:
        --------
        pd.DataFrame
            AUC scores (gene sets x cells)
        """
        
        if rankings is None:
            if self.rankings is None:
                raise ValueError("No rankings available. Run build_rankings first.")
            rankings = self.rankings
            
        print("Calculating AUC scores...")
        
        n_genes = rankings.shape[0]
        
        # Determine max rank threshold
        if isinstance(auc_max_rank, float) and auc_max_rank <= 1.0:
            max_rank = int(n_genes * auc_max_rank)
        else:
            max_rank = int(auc_max_rank)
            
        max_rank = min(max_rank, n_genes)
        
        auc_scores = {}
        
        for gene_set_name, genes in gene_sets.items():
            # Find genes present in rankings
            present_genes = [g for g in genes if g in rankings.index]
            
            if len(present_genes) == 0:
                # No genes found, set AUC to 0
                auc_scores[gene_set_name] = np.zeros(rankings.shape[1])
                continue
                
            # Calculate AUC for each cell
            cell_aucs = []
            
            for cell in rankings.columns:
                cell_rankings = rankings[cell]
                
                # Get ranks of gene set genes
                gene_set_ranks = []
                for gene in present_genes:
                    if gene in cell_rankings.index:
                        rank = cell_rankings[gene]
                        if rank <= max_rank:
                            gene_set_ranks.append(rank)
                
                if len(gene_set_ranks) == 0:
                    auc = 0.0
                else:
                    # Calculate AUC using trapezoidal rule
                    gene_set_ranks = sorted(gene_set_ranks)
                    
                    # Create x and y for AUC calculation
                    x = np.arange(1, len(gene_set_ranks) + 1) / len(present_genes)
                    y = (max_rank - np.array(gene_set_ranks) + 1) / max_rank
                    
                    # Calculate AUC
                    auc = np.trapz(y, x)
                
                cell_aucs.append(auc)
            
            auc_scores[gene_set_name] = np.array(cell_aucs)
        
        # Create DataFrame
        auc_df = pd.DataFrame(auc_scores, index=rankings.columns).T
        
        return auc_df
    
    def _plot_ranking_stats(self, rankings: pd.DataFrame):
        """Plot ranking statistics."""
        try:
            import matplotlib.pyplot as plt
            
            # Plot distribution of rankings
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            
            # Histogram of mean rankings per gene
            mean_rankings = rankings.mean(axis=1)
            axes[0].hist(mean_rankings, bins=50, alpha=0.7)
            axes[0].set_xlabel('Mean Ranking')
            axes[0].set_ylabel('Number of Genes')
            axes[0].set_title('Distribution of Mean Gene Rankings')
            
            # Histogram of number of expressed genes per cell
            expressed_per_cell = (rankings > 0).sum(axis=0)
            axes[1].hist(expressed_per_cell, bins=50, alpha=0.7)
            axes[1].set_xlabel('Number of Expressed Genes')
            axes[1].set_ylabel('Number of Cells')
            axes[1].set_title('Distribution of Expressed Genes per Cell')
            
            plt.tight_layout()
            plt.show()
            
        except ImportError:
            print("Matplotlib not available for plotting")
    
    def get_auc_threshold(
        self,
        auc_scores: pd.DataFrame,
        method: str = "global",
        threshold_quantile: float = 0.01
    ) -> Dict[str, float]:
        """
        Calculate AUC thresholds for significance.
        
        Parameters:
        -----------
        auc_scores : pd.DataFrame
            AUC scores matrix
        method : str
            Method for threshold calculation: "global" or "per_gene_set"
        threshold_quantile : float
            Quantile for threshold calculation
            
        Returns:
        --------
        Dict[str, float]
            Thresholds for each gene set
        """
        
        thresholds = {}
        
        if method == "global":
            # Use global threshold
            global_threshold = np.quantile(auc_scores.values.flatten(), 1 - threshold_quantile)
            for gene_set in auc_scores.index:
                thresholds[gene_set] = global_threshold
                
        elif method == "per_gene_set":
            # Calculate threshold per gene set
            for gene_set in auc_scores.index:
                gene_set_scores = auc_scores.loc[gene_set]
                threshold = np.quantile(gene_set_scores, 1 - threshold_quantile)
                thresholds[gene_set] = threshold
                
        return thresholds