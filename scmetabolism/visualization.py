"""
Visualization functions for metabolism analysis.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, List, Union, Tuple
import warnings


class MetabolismPlotter:
    """
    Visualization class for metabolism analysis results.
    """
    
    def __init__(self):
        self.default_colors = plt.cm.viridis
        
    def dim_plot(
        self,
        embedding: pd.DataFrame,
        metabolism_scores: pd.DataFrame,
        pathway: str,
        embedding_type: str = "umap",
        size: float = 1.0,
        alpha: float = 0.8,
        cmap: str = "viridis",
        figsize: Tuple[int, int] = (8, 6)
    ) -> plt.Figure:
        """
        Create dimension reduction plot colored by metabolism pathway activity.
        
        Parameters:
        -----------
        embedding : pd.DataFrame
            2D embedding coordinates (cells x 2)
        metabolism_scores : pd.DataFrame
            Metabolism scores (pathways x cells)
        pathway : str
            Pathway name to visualize
        embedding_type : str
            Type of embedding ("umap", "tsne", "pca")
        size : float
            Point size
        alpha : float
            Point transparency
        cmap : str
            Colormap name
        figsize : tuple
            Figure size
            
        Returns:
        --------
        plt.Figure
            Matplotlib figure object
        """
        
        if pathway not in metabolism_scores.index:
            raise ValueError(f"Pathway '{pathway}' not found in metabolism scores")
            
        # Get pathway scores
        pathway_scores = metabolism_scores.loc[pathway]
        
        # Align data
        common_cells = embedding.index.intersection(pathway_scores.index)
        embedding_aligned = embedding.loc[common_cells]
        scores_aligned = pathway_scores.loc[common_cells]
        
        # Create plot
        fig, ax = plt.subplots(figsize=figsize)
        
        scatter = ax.scatter(
            embedding_aligned.iloc[:, 0],
            embedding_aligned.iloc[:, 1],
            c=scores_aligned,
            s=size * 20,
            alpha=alpha,
            cmap=cmap,
            edgecolors='none'
        )
        
        # Customize plot
        ax.set_xlabel(f'{embedding_type.upper()} 1')
        ax.set_ylabel(f'{embedding_type.upper()} 2')
        ax.set_title(f'{pathway}')
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Metabolism Score')
        
        # Remove spines
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        
        plt.tight_layout()
        
        print("\nPlease Cite:")
        print("Yingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021.")
        print("https://pubmed.ncbi.nlm.nih.gov/34417225/\n")
        
        return fig
    
    def dot_plot(
        self,
        metabolism_scores: pd.DataFrame,
        metadata: pd.DataFrame,
        pathways: List[str],
        group_by: str,
        normalize: str = "row",
        figsize: Tuple[int, int] = (10, 8),
        cmap: str = "viridis"
    ) -> plt.Figure:
        """
        Create dot plot showing metabolism pathway activity across cell groups.
        
        Parameters:
        -----------
        metabolism_scores : pd.DataFrame
            Metabolism scores (pathways x cells)
        metadata : pd.DataFrame
            Cell metadata
        pathways : List[str]
            List of pathways to plot
        group_by : str
            Column in metadata to group cells by
        normalize : str
            Normalization method: "row", "column", or "none"
        figsize : tuple
            Figure size
        cmap : str
            Colormap name
            
        Returns:
        --------
        plt.Figure
            Matplotlib figure object
        """
        
        # Validate inputs
        missing_pathways = [p for p in pathways if p not in metabolism_scores.index]
        if missing_pathways:
            raise ValueError(f"Pathways not found: {missing_pathways}")
            
        if group_by not in metadata.columns:
            raise ValueError(f"Column '{group_by}' not found in metadata")
        
        # Align data
        common_cells = metabolism_scores.columns.intersection(metadata.index)
        scores_aligned = metabolism_scores.loc[pathways, common_cells]
        metadata_aligned = metadata.loc[common_cells]
        
        # Calculate mean scores per group
        groups = metadata_aligned[group_by].unique()
        mean_scores = []
        
        for group in groups:
            group_cells = metadata_aligned[metadata_aligned[group_by] == group].index
            group_scores = scores_aligned.loc[:, group_cells].mean(axis=1)
            mean_scores.append(group_scores)
        
        # Create matrix for plotting
        plot_matrix = pd.DataFrame(mean_scores, index=groups, columns=pathways).T
        
        # Normalize if requested
        if normalize == "row":
            plot_matrix = plot_matrix.div(plot_matrix.max(axis=1), axis=0)
        elif normalize == "column":
            plot_matrix = plot_matrix.div(plot_matrix.max(axis=0), axis=1)
        
        # Create plot
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create dot plot
        x_pos = np.arange(len(groups))
        y_pos = np.arange(len(pathways))
        
        for i, pathway in enumerate(pathways):
            for j, group in enumerate(groups):
                value = plot_matrix.loc[pathway, group]
                size = abs(value) * 200  # Scale size
                
                ax.scatter(j, i, s=size, c=value, cmap=cmap, 
                          alpha=0.8, edgecolors='black', linewidth=0.5)
        
        # Customize plot
        ax.set_xticks(x_pos)
        ax.set_xticklabels(groups, rotation=45, ha='right')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(pathways)
        ax.set_xlabel(group_by)
        ax.set_ylabel('Metabolic Pathways')
        ax.set_title('Metabolism Pathway Activity')
        
        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap)
        sm.set_array(plot_matrix.values)
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label('Normalized Activity')
        
        plt.tight_layout()
        
        print("\nPlease Cite:")
        print("Yingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021.")
        print("https://pubmed.ncbi.nlm.nih.gov/34417225/\n")
        
        return fig
    
    def box_plot(
        self,
        metabolism_scores: pd.DataFrame,
        metadata: pd.DataFrame,
        pathways: List[str],
        group_by: str,
        ncols: int = 2,
        figsize: Tuple[int, int] = (12, 8)
    ) -> plt.Figure:
        """
        Create box plots showing metabolism pathway activity distributions.
        
        Parameters:
        -----------
        metabolism_scores : pd.DataFrame
            Metabolism scores (pathways x cells)
        metadata : pd.DataFrame
            Cell metadata
        pathways : List[str]
            List of pathways to plot
        group_by : str
            Column in metadata to group cells by
        ncols : int
            Number of columns in subplot grid
        figsize : tuple
            Figure size
            
        Returns:
        --------
        plt.Figure
            Matplotlib figure object
        """
        
        # Validate inputs
        missing_pathways = [p for p in pathways if p not in metabolism_scores.index]
        if missing_pathways:
            raise ValueError(f"Pathways not found: {missing_pathways}")
            
        if group_by not in metadata.columns:
            raise ValueError(f"Column '{group_by}' not found in metadata")
        
        # Align data
        common_cells = metabolism_scores.columns.intersection(metadata.index)
        scores_aligned = metabolism_scores.loc[pathways, common_cells]
        metadata_aligned = metadata.loc[common_cells]
        
        # Calculate subplot layout
        nrows = int(np.ceil(len(pathways) / ncols))
        
        # Create subplots
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        if nrows == 1:
            axes = [axes] if ncols == 1 else axes
        else:
            axes = axes.flatten()
        
        # Create box plots
        for i, pathway in enumerate(pathways):
            ax = axes[i]
            
            # Prepare data for plotting
            plot_data = []
            plot_labels = []
            
            for group in metadata_aligned[group_by].unique():
                group_cells = metadata_aligned[metadata_aligned[group_by] == group].index
                group_scores = scores_aligned.loc[pathway, group_cells]
                plot_data.append(group_scores.values)
                plot_labels.append(group)
            
            # Create box plot
            bp = ax.boxplot(plot_data, labels=plot_labels, patch_artist=True)
            
            # Color boxes
            colors = plt.cm.Set3(np.linspace(0, 1, len(plot_labels)))
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            
            ax.set_title(pathway)
            ax.set_xlabel(group_by)
            ax.set_ylabel('Metabolism Score')
            
            # Rotate x-axis labels if needed
            if len(plot_labels) > 3:
                ax.tick_params(axis='x', rotation=45)
        
        # Hide unused subplots
        for i in range(len(pathways), len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        
        print("\nPlease Cite:")
        print("Yingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021.")
        print("https://pubmed.ncbi.nlm.nih.gov/34417225/\n")
        
        return fig
    
    def heatmap(
        self,
        metabolism_scores: pd.DataFrame,
        pathways: Optional[List[str]] = None,
        cells: Optional[List[str]] = None,
        cluster_rows: bool = True,
        cluster_cols: bool = True,
        figsize: Tuple[int, int] = (12, 8),
        cmap: str = "viridis"
    ) -> plt.Figure:
        """
        Create heatmap of metabolism scores.
        
        Parameters:
        -----------
        metabolism_scores : pd.DataFrame
            Metabolism scores (pathways x cells)
        pathways : List[str], optional
            Subset of pathways to plot
        cells : List[str], optional
            Subset of cells to plot
        cluster_rows : bool
            Whether to cluster rows
        cluster_cols : bool
            Whether to cluster columns
        figsize : tuple
            Figure size
        cmap : str
            Colormap name
            
        Returns:
        --------
        plt.Figure
            Matplotlib figure object
        """
        
        # Subset data if requested
        plot_data = metabolism_scores.copy()
        
        if pathways is not None:
            missing_pathways = [p for p in pathways if p not in plot_data.index]
            if missing_pathways:
                warnings.warn(f"Pathways not found: {missing_pathways}")
            pathways = [p for p in pathways if p in plot_data.index]
            plot_data = plot_data.loc[pathways]
        
        if cells is not None:
            missing_cells = [c for c in cells if c not in plot_data.columns]
            if missing_cells:
                warnings.warn(f"Cells not found: {missing_cells}")
            cells = [c for c in cells if c in plot_data.columns]
            plot_data = plot_data.loc[:, cells]
        
        # Create heatmap
        fig, ax = plt.subplots(figsize=figsize)
        
        sns.heatmap(
            plot_data,
            cmap=cmap,
            center=0,
            robust=True,
            cbar_kws={'label': 'Metabolism Score'},
            ax=ax
        )
        
        ax.set_title('Metabolism Pathway Activity Heatmap')
        ax.set_xlabel('Cells')
        ax.set_ylabel('Pathways')
        
        plt.tight_layout()
        
        return fig