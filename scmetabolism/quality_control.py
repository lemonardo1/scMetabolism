"""
Quality control and data validation for scMetabolism.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional, Union
import warnings
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


class DataValidator:
    """Validate input data for scMetabolism analysis."""
    
    def __init__(self):
        self.validation_results = {}
    
    def validate_count_matrix(
        self,
        count_matrix: pd.DataFrame,
        min_genes_per_cell: int = 200,
        min_cells_per_gene: int = 3,
        max_mito_percent: float = 20.0
    ) -> Dict[str, Union[bool, str, int, float]]:
        """
        Validate single-cell count matrix.
        
        Parameters:
        -----------
        count_matrix : pd.DataFrame
            Gene expression count matrix (genes x cells)
        min_genes_per_cell : int
            Minimum genes per cell threshold
        min_cells_per_gene : int
            Minimum cells per gene threshold
        max_mito_percent : float
            Maximum mitochondrial gene percentage
            
        Returns:
        --------
        Dict
            Validation results
        """
        
        results = {
            'is_valid': True,
            'warnings': [],
            'errors': [],
            'stats': {}
        }
        
        # Basic shape validation
        n_genes, n_cells = count_matrix.shape
        results['stats']['n_genes'] = n_genes
        results['stats']['n_cells'] = n_cells
        
        if n_genes == 0 or n_cells == 0:
            results['errors'].append("Empty count matrix")
            results['is_valid'] = False
            return results
        
        # Check for negative values
        if (count_matrix < 0).any().any():
            results['errors'].append("Negative values found in count matrix")
            results['is_valid'] = False
        
        # Check data type
        if not np.issubdtype(count_matrix.dtypes.iloc[0], np.number):
            results['errors'].append("Non-numeric data found")
            results['is_valid'] = False
        
        # Gene statistics
        genes_per_cell = (count_matrix > 0).sum(axis=0)
        cells_per_gene = (count_matrix > 0).sum(axis=1)
        
        results['stats']['mean_genes_per_cell'] = genes_per_cell.mean()
        results['stats']['median_genes_per_cell'] = genes_per_cell.median()
        results['stats']['mean_cells_per_gene'] = cells_per_gene.mean()
        results['stats']['median_cells_per_gene'] = cells_per_gene.median()
        
        # Low quality cells
        low_gene_cells = (genes_per_cell < min_genes_per_cell).sum()
        if low_gene_cells > 0:
            results['warnings'].append(
                f"{low_gene_cells} cells have < {min_genes_per_cell} genes"
            )
        
        # Low expression genes
        low_expr_genes = (cells_per_gene < min_cells_per_gene).sum()
        if low_expr_genes > 0:
            results['warnings'].append(
                f"{low_expr_genes} genes expressed in < {min_cells_per_gene} cells"
            )
        
        # Mitochondrial gene percentage
        mito_genes = [gene for gene in count_matrix.index 
                     if gene.startswith(('MT-', 'mt-', 'Mt-'))]
        
        if len(mito_genes) > 0:
            mito_counts = count_matrix.loc[mito_genes].sum(axis=0)
            total_counts = count_matrix.sum(axis=0)
            mito_percent = (mito_counts / total_counts * 100).fillna(0)
            
            results['stats']['n_mito_genes'] = len(mito_genes)
            results['stats']['mean_mito_percent'] = mito_percent.mean()
            results['stats']['median_mito_percent'] = mito_percent.median()
            
            high_mito_cells = (mito_percent > max_mito_percent).sum()
            if high_mito_cells > 0:
                results['warnings'].append(
                    f"{high_mito_cells} cells have > {max_mito_percent}% mitochondrial genes"
                )
        
        # Sparsity
        sparsity = (count_matrix == 0).sum().sum() / (n_genes * n_cells)
        results['stats']['sparsity'] = sparsity
        
        if sparsity > 0.95:
            results['warnings'].append(f"Very sparse matrix ({sparsity:.1%} zeros)")
        
        # Total UMI distribution
        total_umi = count_matrix.sum(axis=0)
        results['stats']['mean_umi_per_cell'] = total_umi.mean()
        results['stats']['median_umi_per_cell'] = total_umi.median()
        results['stats']['std_umi_per_cell'] = total_umi.std()
        
        # Check for outliers
        q1, q3 = total_umi.quantile([0.25, 0.75])
        iqr = q3 - q1
        outlier_threshold_low = q1 - 1.5 * iqr
        outlier_threshold_high = q3 + 1.5 * iqr
        
        outlier_cells = ((total_umi < outlier_threshold_low) | 
                        (total_umi > outlier_threshold_high)).sum()
        
        if outlier_cells > 0:
            results['warnings'].append(f"{outlier_cells} potential outlier cells detected")
        
        self.validation_results = results
        return results
    
    def validate_gene_sets(
        self,
        gene_sets: Dict[str, List[str]],
        count_matrix: pd.DataFrame
    ) -> Dict[str, Union[bool, List, Dict]]:
        """
        Validate gene sets against count matrix.
        
        Parameters:
        -----------
        gene_sets : Dict[str, List[str]]
            Gene sets dictionary
        count_matrix : pd.DataFrame
            Count matrix
            
        Returns:
        --------
        Dict
            Validation results
        """
        
        results = {
            'is_valid': True,
            'warnings': [],
            'pathway_stats': {},
            'overall_stats': {}
        }
        
        available_genes = set(count_matrix.index)
        total_pathways = len(gene_sets)
        valid_pathways = 0
        total_genes_in_sets = 0
        total_matched_genes = 0
        
        for pathway_name, genes in gene_sets.items():
            pathway_genes = set(genes)
            matched_genes = pathway_genes.intersection(available_genes)
            
            total_genes_in_sets += len(pathway_genes)
            total_matched_genes += len(matched_genes)
            
            match_rate = len(matched_genes) / len(pathway_genes) if pathway_genes else 0
            
            results['pathway_stats'][pathway_name] = {
                'total_genes': len(pathway_genes),
                'matched_genes': len(matched_genes),
                'match_rate': match_rate
            }
            
            if len(matched_genes) >= 2:  # Minimum genes for analysis
                valid_pathways += 1
            elif len(matched_genes) == 0:
                results['warnings'].append(f"No genes found for pathway: {pathway_name}")
            else:
                results['warnings'].append(f"Only {len(matched_genes)} gene found for pathway: {pathway_name}")
        
        # Overall statistics
        overall_match_rate = total_matched_genes / total_genes_in_sets if total_genes_in_sets > 0 else 0
        
        results['overall_stats'] = {
            'total_pathways': total_pathways,
            'valid_pathways': valid_pathways,
            'total_genes_in_sets': total_genes_in_sets,
            'total_matched_genes': total_matched_genes,
            'overall_match_rate': overall_match_rate
        }
        
        if valid_pathways == 0:
            results['is_valid'] = False
            results['warnings'].append("No valid pathways found")
        elif valid_pathways < total_pathways * 0.5:
            results['warnings'].append(f"Only {valid_pathways}/{total_pathways} pathways are valid")
        
        return results
    
    def plot_validation_summary(self, figsize: Tuple[int, int] = (15, 10)):
        """Plot validation summary."""
        
        if not self.validation_results:
            print("No validation results available. Run validate_count_matrix first.")
            return
        
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        fig.suptitle('Data Quality Control Summary', fontsize=16)
        
        stats = self.validation_results['stats']
        
        # 1. Genes per cell distribution
        if 'mean_genes_per_cell' in stats:
            ax = axes[0, 0]
            # This would need actual data to plot histogram
            ax.text(0.5, 0.5, f"Mean genes/cell:\n{stats['mean_genes_per_cell']:.0f}", 
                   ha='center', va='center', transform=ax.transAxes, fontsize=12)
            ax.set_title('Genes per Cell')
        
        # 2. UMI per cell distribution  
        if 'mean_umi_per_cell' in stats:
            ax = axes[0, 1]
            ax.text(0.5, 0.5, f"Mean UMI/cell:\n{stats['mean_umi_per_cell']:.0f}", 
                   ha='center', va='center', transform=ax.transAxes, fontsize=12)
            ax.set_title('UMI per Cell')
        
        # 3. Mitochondrial percentage
        if 'mean_mito_percent' in stats:
            ax = axes[0, 2]
            ax.text(0.5, 0.5, f"Mean mito %:\n{stats['mean_mito_percent']:.1f}%", 
                   ha='center', va='center', transform=ax.transAxes, fontsize=12)
            ax.set_title('Mitochondrial %')
        
        # 4. Sparsity
        ax = axes[1, 0]
        ax.text(0.5, 0.5, f"Sparsity:\n{stats['sparsity']:.1%}", 
               ha='center', va='center', transform=ax.transAxes, fontsize=12)
        ax.set_title('Matrix Sparsity')
        
        # 5. Data dimensions
        ax = axes[1, 1]
        ax.text(0.5, 0.5, f"Dimensions:\n{stats['n_genes']} genes\n{stats['n_cells']} cells", 
               ha='center', va='center', transform=ax.transAxes, fontsize=12)
        ax.set_title('Data Dimensions')
        
        # 6. Warnings/Errors
        ax = axes[1, 2]
        n_warnings = len(self.validation_results['warnings'])
        n_errors = len(self.validation_results['errors'])
        
        status_text = f"Warnings: {n_warnings}\nErrors: {n_errors}"
        color = 'green' if n_errors == 0 else 'red'
        
        ax.text(0.5, 0.5, status_text, ha='center', va='center', 
               transform=ax.transAxes, fontsize=12, color=color)
        ax.set_title('Validation Status')
        
        # Remove axes
        for ax in axes.flat:
            ax.set_xticks([])
            ax.set_yticks([])
        
        plt.tight_layout()
        return fig


class QualityMetrics:
    """Calculate quality metrics for metabolism analysis."""
    
    @staticmethod
    def calculate_pathway_coverage(
        gene_sets: Dict[str, List[str]],
        count_matrix: pd.DataFrame
    ) -> pd.DataFrame:
        """Calculate pathway coverage statistics."""
        
        coverage_stats = []
        available_genes = set(count_matrix.index)
        
        for pathway_name, genes in gene_sets.items():
            pathway_genes = set(genes)
            matched_genes = pathway_genes.intersection(available_genes)
            
            coverage_stats.append({
                'pathway': pathway_name,
                'total_genes': len(pathway_genes),
                'matched_genes': len(matched_genes),
                'coverage_rate': len(matched_genes) / len(pathway_genes) if pathway_genes else 0,
                'missing_genes': len(pathway_genes) - len(matched_genes)
            })
        
        return pd.DataFrame(coverage_stats)
    
    @staticmethod
    def calculate_score_quality(
        metabolism_scores: pd.DataFrame
    ) -> Dict[str, float]:
        """Calculate quality metrics for metabolism scores."""
        
        metrics = {}
        
        # Score distribution metrics
        metrics['mean_score'] = metabolism_scores.values.mean()
        metrics['std_score'] = metabolism_scores.values.std()
        metrics['min_score'] = metabolism_scores.values.min()
        metrics['max_score'] = metabolism_scores.values.max()
        
        # Pathway-level metrics
        pathway_means = metabolism_scores.mean(axis=1)
        pathway_stds = metabolism_scores.std(axis=1)
        
        metrics['pathway_mean_range'] = pathway_means.max() - pathway_means.min()
        metrics['pathway_std_range'] = pathway_stds.max() - pathway_stds.min()
        
        # Cell-level metrics
        cell_means = metabolism_scores.mean(axis=0)
        cell_stds = metabolism_scores.std(axis=0)
        
        metrics['cell_mean_range'] = cell_means.max() - cell_means.min()
        metrics['cell_std_range'] = cell_stds.max() - cell_stds.min()
        
        # Correlation structure
        pathway_corr = metabolism_scores.T.corr()
        metrics['mean_pathway_correlation'] = pathway_corr.values[np.triu_indices_from(pathway_corr.values, k=1)].mean()
        
        return metrics
    
    @staticmethod
    def detect_outliers(
        metabolism_scores: pd.DataFrame,
        method: str = "iqr",
        threshold: float = 1.5
    ) -> Dict[str, List[str]]:
        """Detect outlier cells or pathways."""
        
        outliers = {'cells': [], 'pathways': []}
        
        if method == "iqr":
            # Cell outliers
            cell_means = metabolism_scores.mean(axis=0)
            q1, q3 = cell_means.quantile([0.25, 0.75])
            iqr = q3 - q1
            lower_bound = q1 - threshold * iqr
            upper_bound = q3 + threshold * iqr
            
            outlier_cells = cell_means[(cell_means < lower_bound) | (cell_means > upper_bound)]
            outliers['cells'] = outlier_cells.index.tolist()
            
            # Pathway outliers
            pathway_means = metabolism_scores.mean(axis=1)
            q1, q3 = pathway_means.quantile([0.25, 0.75])
            iqr = q3 - q1
            lower_bound = q1 - threshold * iqr
            upper_bound = q3 + threshold * iqr
            
            outlier_pathways = pathway_means[(pathway_means < lower_bound) | (pathway_means > upper_bound)]
            outliers['pathways'] = outlier_pathways.index.tolist()
        
        elif method == "zscore":
            # Z-score based outlier detection
            cell_means = metabolism_scores.mean(axis=0)
            cell_zscores = np.abs(stats.zscore(cell_means))
            outliers['cells'] = cell_means[cell_zscores > threshold].index.tolist()
            
            pathway_means = metabolism_scores.mean(axis=1)
            pathway_zscores = np.abs(stats.zscore(pathway_means))
            outliers['pathways'] = pathway_means[pathway_zscores > threshold].index.tolist()
        
        return outliers