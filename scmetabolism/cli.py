"""
Command Line Interface for scMetabolism.
"""

import argparse
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import json
import yaml

from .core import ScMetabolism
from .visualization import MetabolismPlotter
from .quality_control import DataValidator
from .config import get_config


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="scMetabolism: Single-cell metabolism analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis
  scmetabolism analyze --input data.csv --output results/ --method aucell
  
  # GO analysis
  scmetabolism analyze --input data.csv --gene-sets GO_metabolism --output results/
  
  # Quality control
  scmetabolism qc --input data.csv --output qc_report.html
  
  # Configuration
  scmetabolism config --show
  scmetabolism config --set analysis.default_method=ssgsea
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Analyze command
    analyze_parser = subparsers.add_parser('analyze', help='Run metabolism analysis')
    analyze_parser.add_argument('--input', '-i', required=True, 
                               help='Input count matrix file (CSV/TSV/H5AD)')
    analyze_parser.add_argument('--output', '-o', required=True,
                               help='Output directory')
    analyze_parser.add_argument('--method', default='aucell',
                               choices=['aucell', 'ssgsea', 'gsva'],
                               help='Analysis method')
    analyze_parser.add_argument('--gene-sets', default='KEGG',
                               help='Gene sets to use (KEGG, REACTOME, GO_metabolism, etc.)')
    analyze_parser.add_argument('--imputation', action='store_true',
                               help='Enable ALRA imputation')
    analyze_parser.add_argument('--n-cores', type=int, default=2,
                               help='Number of cores to use')
    analyze_parser.add_argument('--metadata', 
                               help='Cell metadata file (CSV/TSV)')
    analyze_parser.add_argument('--embedding',
                               help='2D embedding coordinates file (CSV/TSV)')
    
    # Quality control command
    qc_parser = subparsers.add_parser('qc', help='Run quality control')
    qc_parser.add_argument('--input', '-i', required=True,
                          help='Input count matrix file')
    qc_parser.add_argument('--output', '-o', required=True,
                          help='Output report file (HTML/PDF)')
    qc_parser.add_argument('--min-genes', type=int, default=200,
                          help='Minimum genes per cell')
    qc_parser.add_argument('--min-cells', type=int, default=3,
                          help='Minimum cells per gene')
    qc_parser.add_argument('--max-mito', type=float, default=20.0,
                          help='Maximum mitochondrial percentage')
    
    # Configuration command
    config_parser = subparsers.add_parser('config', help='Manage configuration')
    config_parser.add_argument('--show', action='store_true',
                              help='Show current configuration')
    config_parser.add_argument('--set', 
                              help='Set configuration value (key=value)')
    config_parser.add_argument('--reset', action='store_true',
                              help='Reset to default configuration')
    config_parser.add_argument('--template', action='store_true',
                              help='Create configuration template')
    
    # Version command
    version_parser = subparsers.add_parser('version', help='Show version')
    
    args = parser.parse_args()
    
    if args.command == 'analyze':
        run_analysis(args)
    elif args.command == 'qc':
        run_quality_control(args)
    elif args.command == 'config':
        manage_config(args)
    elif args.command == 'version':
        show_version()
    else:
        parser.print_help()


def run_analysis(args):
    """Run metabolism analysis."""
    print(f"🔬 Starting scMetabolism analysis...")
    print(f"   Input: {args.input}")
    print(f"   Method: {args.method}")
    print(f"   Gene sets: {args.gene_sets}")
    
    # Load data
    print("📊 Loading data...")
    count_matrix = load_data(args.input)
    
    # Load metadata if provided
    metadata = None
    if args.metadata:
        metadata = pd.read_csv(args.metadata, index_col=0)
        print(f"   Loaded metadata: {metadata.shape}")
    
    # Load embedding if provided
    embedding = None
    if args.embedding:
        embedding = pd.read_csv(args.embedding, index_col=0)
        print(f"   Loaded embedding: {embedding.shape}")
    
    # Run analysis
    print("🧬 Computing metabolism scores...")
    sc_metab = ScMetabolism()
    
    metabolism_scores = sc_metab.compute_metabolism(
        count_matrix=count_matrix,
        method=args.method,
        imputation=args.imputation,
        n_cores=args.n_cores,
        metabolism_type=args.gene_sets
    )
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save results
    print("💾 Saving results...")
    metabolism_scores.to_csv(output_dir / "metabolism_scores.csv")
    
    # Create visualizations
    if metadata is not None or embedding is not None:
        print("📈 Creating visualizations...")
        plotter = MetabolismPlotter()
        
        # Get top pathways
        mean_scores = metabolism_scores.mean(axis=1).sort_values(ascending=False)
        top_pathways = mean_scores.head(6).index.tolist()
        
        if embedding is not None and len(top_pathways) > 0:
            # Dimension plots
            for pathway in top_pathways[:3]:
                try:
                    fig = plotter.dim_plot(
                        embedding=embedding,
                        metabolism_scores=metabolism_scores,
                        pathway=pathway,
                        embedding_type="embedding"
                    )
                    fig.savefig(output_dir / f"dimplot_{pathway.replace('/', '_')}.png", 
                               dpi=300, bbox_inches='tight')
                    plt.close(fig)
                except Exception as e:
                    print(f"   Warning: Could not create plot for {pathway}: {e}")
        
        if metadata is not None and len(top_pathways) > 0:
            # Box plots
            try:
                fig = plotter.box_plot(
                    metabolism_scores=metabolism_scores,
                    metadata=metadata,
                    pathways=top_pathways,
                    group_by=metadata.columns[0],  # Use first metadata column
                    ncols=2
                )
                fig.savefig(output_dir / "boxplots.png", dpi=300, bbox_inches='tight')
                plt.close(fig)
            except Exception as e:
                print(f"   Warning: Could not create box plots: {e}")
        
        # Heatmap
        try:
            fig = plotter.heatmap(
                metabolism_scores=metabolism_scores,
                pathways=top_pathways,
                cluster_rows=True,
                cluster_cols=False
            )
            fig.savefig(output_dir / "heatmap.png", dpi=300, bbox_inches='tight')
            plt.close(fig)
        except Exception as e:
            print(f"   Warning: Could not create heatmap: {e}")
    
    # Save summary
    summary = {
        'analysis_method': args.method,
        'gene_sets': args.gene_sets,
        'imputation': args.imputation,
        'n_pathways': metabolism_scores.shape[0],
        'n_cells': metabolism_scores.shape[1],
        'top_pathways': mean_scores.head(10).to_dict()
    }
    
    with open(output_dir / "analysis_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"✅ Analysis complete! Results saved to: {output_dir}")
    print(f"   - Metabolism scores: metabolism_scores.csv")
    print(f"   - Analysis summary: analysis_summary.json")
    if metadata is not None or embedding is not None:
        print(f"   - Visualizations: *.png files")


def run_quality_control(args):
    """Run quality control analysis."""
    print(f"🔍 Running quality control...")
    print(f"   Input: {args.input}")
    
    # Load data
    count_matrix = load_data(args.input)
    
    # Run validation
    validator = DataValidator()
    results = validator.validate_count_matrix(
        count_matrix,
        min_genes_per_cell=args.min_genes,
        min_cells_per_gene=args.min_cells,
        max_mito_percent=args.max_mito
    )
    
    # Create report
    output_file = Path(args.output)
    
    if output_file.suffix.lower() == '.html':
        create_html_report(results, output_file)
    else:
        # Save as JSON
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
    
    print(f"✅ Quality control complete! Report saved to: {output_file}")


def manage_config(args):
    """Manage configuration."""
    config = get_config()
    
    if args.show:
        print("📋 Current configuration:")
        print(yaml.dump(config.config, default_flow_style=False))
    
    elif args.set:
        key, value = args.set.split('=', 1)
        # Try to parse value as appropriate type
        try:
            if value.lower() in ['true', 'false']:
                value = value.lower() == 'true'
            elif value.isdigit():
                value = int(value)
            elif '.' in value and value.replace('.', '').isdigit():
                value = float(value)
        except:
            pass  # Keep as string
        
        config.set(key, value)
        config.save_config()
        print(f"✅ Set {key} = {value}")
    
    elif args.reset:
        config.reset_to_defaults()
        config.save_config()
        print("✅ Configuration reset to defaults")
    
    elif args.template:
        config.create_user_config_template()
        print("✅ Configuration template created")


def show_version():
    """Show version information."""
    from . import __version__
    print(f"scMetabolism version {__version__}")


def load_data(filepath):
    """Load data from various formats."""
    filepath = Path(filepath)
    
    if filepath.suffix.lower() == '.h5ad':
        try:
            import scanpy as sc
            adata = sc.read_h5ad(filepath)
            return pd.DataFrame(
                adata.X.T if hasattr(adata.X, 'todense') else adata.X.T,
                index=adata.var_names,
                columns=adata.obs_names
            )
        except ImportError:
            raise ImportError("scanpy required for H5AD files. Install with: pip install scanpy")
    
    elif filepath.suffix.lower() in ['.csv', '.tsv', '.txt']:
        sep = ',' if filepath.suffix.lower() == '.csv' else '\t'
        return pd.read_csv(filepath, index_col=0, sep=sep)
    
    else:
        raise ValueError(f"Unsupported file format: {filepath.suffix}")


def create_html_report(results, output_file):
    """Create HTML quality control report."""
    html_template = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>scMetabolism Quality Control Report</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            .header { background-color: #f0f0f0; padding: 20px; border-radius: 5px; }
            .section { margin: 20px 0; }
            .stats { background-color: #f9f9f9; padding: 15px; border-radius: 5px; }
            .warning { color: orange; }
            .error { color: red; }
            .success { color: green; }
        </style>
    </head>
    <body>
        <div class="header">
            <h1>scMetabolism Quality Control Report</h1>
            <p>Generated on: {timestamp}</p>
        </div>
        
        <div class="section">
            <h2>Validation Status</h2>
            <p class="{status_class}">Status: {status}</p>
        </div>
        
        <div class="section">
            <h2>Data Statistics</h2>
            <div class="stats">
                {stats_html}
            </div>
        </div>
        
        <div class="section">
            <h2>Warnings</h2>
            {warnings_html}
        </div>
        
        <div class="section">
            <h2>Errors</h2>
            {errors_html}
        </div>
    </body>
    </html>
    """
    
    from datetime import datetime
    
    # Format statistics
    stats_html = ""
    for key, value in results['stats'].items():
        if isinstance(value, float):
            value = f"{value:.2f}"
        stats_html += f"<p><strong>{key.replace('_', ' ').title()}:</strong> {value}</p>\n"
    
    # Format warnings
    warnings_html = ""
    if results['warnings']:
        for warning in results['warnings']:
            warnings_html += f'<p class="warning">⚠️ {warning}</p>\n'
    else:
        warnings_html = '<p class="success">✅ No warnings</p>'
    
    # Format errors
    errors_html = ""
    if results['errors']:
        for error in results['errors']:
            errors_html += f'<p class="error">❌ {error}</p>\n'
    else:
        errors_html = '<p class="success">✅ No errors</p>'
    
    # Fill template
    html_content = html_template.format(
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        status="PASS" if results['is_valid'] else "FAIL",
        status_class="success" if results['is_valid'] else "error",
        stats_html=stats_html,
        warnings_html=warnings_html,
        errors_html=errors_html
    )
    
    with open(output_file, 'w') as f:
        f.write(html_content)


if __name__ == '__main__':
    main()