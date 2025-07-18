"""
Visualization functions for metabolism analysis.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, List, Union, Tuple
import warnings
try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    
try:
    import scanpy as sc
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False


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
 
   def interactive_dim_plot(
        self,
        embedding: pd.DataFrame,
        metabolism_scores: pd.DataFrame,
        pathway: str,
        metadata: Optional[pd.DataFrame] = None,
        color_by: Optional[str] = None,
        embedding_type: str = "umap",
        size: int = 5,
        width: int = 800,
        height: int = 600
    ):
        """
        Create interactive dimension reduction plot using Plotly.
        
        Parameters:
        -----------
        embedding : pd.DataFrame
            2D embedding coordinates
        metabolism_scores : pd.DataFrame
            Metabolism scores
        pathway : str
            Pathway to visualize
        metadata : pd.DataFrame, optional
            Cell metadata for additional coloring
        color_by : str, optional
            Column in metadata to color by
        embedding_type : str
            Type of embedding
        size : int
            Point size
        width, height : int
            Plot dimensions
            
        Returns:
        --------
        plotly.graph_objects.Figure or None
            Interactive plot (if plotly available)
        """
        
        if not PLOTLY_AVAILABLE:
            warnings.warn("Plotly not available. Install with: pip install plotly")
            return None
            
        if pathway not in metabolism_scores.index:
            raise ValueError(f"Pathway '{pathway}' not found in metabolism scores")
        
        # Align data
        common_cells = embedding.index.intersection(metabolism_scores.columns)
        embedding_aligned = embedding.loc[common_cells]
        pathway_scores = metabolism_scores.loc[pathway, common_cells]
        
        # Prepare data for plotting
        plot_data = pd.DataFrame({
            f'{embedding_type.upper()}_1': embedding_aligned.iloc[:, 0],
            f'{embedding_type.upper()}_2': embedding_aligned.iloc[:, 1],
            'Metabolism_Score': pathway_scores,
            'Cell_ID': common_cells
        })
        
        # Add metadata if provided
        if metadata is not None and color_by is not None:
            if color_by in metadata.columns:
                metadata_aligned = metadata.loc[common_cells]
                plot_data[color_by] = metadata_aligned[color_by]
                color_column = color_by
            else:
                warnings.warn(f"Column '{color_by}' not found in metadata")
                color_column = 'Metabolism_Score'
        else:
            color_column = 'Metabolism_Score'
        
        # Create interactive plot
        fig = px.scatter(
            plot_data,
            x=f'{embedding_type.upper()}_1',
            y=f'{embedding_type.upper()}_2',
            color=color_column,
            hover_data=['Cell_ID', 'Metabolism_Score'],
            title=f'{pathway} - {embedding_type.upper()} Plot',
            width=width,
            height=height
        )
        
        fig.update_traces(marker=dict(size=size))
        fig.update_layout(
            xaxis_title=f'{embedding_type.upper()} 1',
            yaxis_title=f'{embedding_type.upper()} 2'
        )
        
        return fig
    
    def interactive_heatmap(
        self,
        metabolism_scores: pd.DataFrame,
        pathways: Optional[List[str]] = None,
        cells: Optional[List[str]] = None,
        metadata: Optional[pd.DataFrame] = None,
        width: int = 1000,
        height: int = 800
    ):
        """
        Create interactive heatmap using Plotly.
        
        Parameters:
        -----------
        metabolism_scores : pd.DataFrame
            Metabolism scores
        pathways : List[str], optional
            Subset of pathways
        cells : List[str], optional
            Subset of cells
        metadata : pd.DataFrame, optional
            Cell metadata for annotations
        width, height : int
            Plot dimensions
            
        Returns:
        --------
        plotly.graph_objects.Figure or None
            Interactive heatmap
        """
        
        if not PLOTLY_AVAILABLE:
            warnings.warn("Plotly not available. Install with: pip install plotly")
            return None
        
        # Subset data
        plot_data = metabolism_scores.copy()
        
        if pathways is not None:
            pathways = [p for p in pathways if p in plot_data.index]
            plot_data = plot_data.loc[pathways]
        
        if cells is not None:
            cells = [c for c in cells if c in plot_data.columns]
            plot_data = plot_data.loc[:, cells]
        
        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=plot_data.values,
            x=plot_data.columns,
            y=plot_data.index,
            colorscale='Viridis',
            hovertemplate='Pathway: %{y}<br>Cell: %{x}<br>Score: %{z}<extra></extra>'
        ))
        
        fig.update_layout(
            title='Interactive Metabolism Heatmap',
            xaxis_title='Cells',
            yaxis_title='Pathways',
            width=width,
            height=height
        )
        
        return fig
    
    def pathway_network_plot(
        self,
        metabolism_scores: pd.DataFrame,
        correlation_threshold: float = 0.5,
        layout: str = "spring",
        node_size_col: Optional[str] = None,
        width: int = 800,
        height: int = 600
    ):
        """
        Create pathway correlation network plot.
        
        Parameters:
        -----------
        metabolism_scores : pd.DataFrame
            Metabolism scores
        correlation_threshold : float
            Minimum correlation to show edge
        layout : str
            Network layout algorithm
        node_size_col : str, optional
            Column to determine node size
        width, height : int
            Plot dimensions
            
        Returns:
        --------
        plotly.graph_objects.Figure or None
            Network plot
        """
        
        if not PLOTLY_AVAILABLE:
            warnings.warn("Plotly not available. Install with: pip install plotly")
            return None
        
        # Calculate pathway correlations
        corr_matrix = metabolism_scores.T.corr()
        
        # Create network edges
        edges = []
        edge_weights = []
        
        for i, pathway1 in enumerate(corr_matrix.index):
            for j, pathway2 in enumerate(corr_matrix.columns):
                if i < j:  # Avoid duplicates
                    corr_val = corr_matrix.loc[pathway1, pathway2]
                    if abs(corr_val) >= correlation_threshold:
                        edges.append((pathway1, pathway2))
                        edge_weights.append(corr_val)
        
        if len(edges) == 0:
            warnings.warn("No correlations above threshold found")
            return None
        
        # Create network layout (simplified)
        import networkx as nx
        
        G = nx.Graph()
        G.add_weighted_edges_from([(e[0], e[1], w) for e, w in zip(edges, edge_weights)])
        
        if layout == "spring":
            pos = nx.spring_layout(G)
        elif layout == "circular":
            pos = nx.circular_layout(G)
        else:
            pos = nx.random_layout(G)
        
        # Create plotly network
        edge_x = []
        edge_y = []
        
        for edge in G.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
        
        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color='#888'),
            hoverinfo='none',
            mode='lines'
        )
        
        node_x = []
        node_y = []
        node_text = []
        
        for node in G.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            node_text.append(node)
        
        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers+text',
            hoverinfo='text',
            text=node_text,
            textposition="middle center",
            marker=dict(
                size=10,
                color='lightblue',
                line=dict(width=2, color='black')
            )
        )
        
        fig = go.Figure(data=[edge_trace, node_trace],
                       layout=go.Layout(
                           title='Pathway Correlation Network',
                           titlefont_size=16,
                           showlegend=False,
                           hovermode='closest',
                           margin=dict(b=20,l=5,r=5,t=40),
                           annotations=[ dict(
                               text=f"Correlation threshold: {correlation_threshold}",
                               showarrow=False,
                               xref="paper", yref="paper",
                               x=0.005, y=-0.002,
                               xanchor="left", yanchor="bottom",
                               font=dict(size=12)
                           )],
                           xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                           yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                           width=width,
                           height=height
                       ))
        
        return fig
    
    def scanpy_integration_plot(
        self,
        adata,
        pathway: str,
        basis: str = "umap",
        layer: Optional[str] = None,
        **kwargs
    ):
        """
        Create scanpy-style plots for metabolism scores.
        
        Parameters:
        -----------
        adata : AnnData
            Annotated data object with metabolism scores
        pathway : str
            Pathway to visualize
        basis : str
            Embedding basis to use
        layer : str, optional
            Layer to use
        **kwargs
            Additional arguments for scanpy plotting
            
        Returns:
        --------
        matplotlib figure or None
        """
        
        if not SCANPY_AVAILABLE:
            warnings.warn("Scanpy not available. Install with: pip install scanpy")
            return None
        
        if 'metabolism' not in adata.obsm:
            raise ValueError("No metabolism scores found in adata.obsm['metabolism']")
        
        # Add pathway score to obs for plotting
        metabolism_df = pd.DataFrame(
            adata.obsm['metabolism'],
            index=adata.obs_names,
            columns=adata.uns['metabolism_pathways']
        )
        
        if pathway not in metabolism_df.columns:
            raise ValueError(f"Pathway '{pathway}' not found in metabolism scores")
        
        adata.obs[f'metabolism_{pathway}'] = metabolism_df[pathway]
        
        # Create scanpy plot
        fig = sc.pl.embedding(
            adata,
            basis=basis,
            color=f'metabolism_{pathway}',
            title=f'{pathway}',
            return_fig=True,
            **kwargs
        )
        
        return fig