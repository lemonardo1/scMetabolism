"""
Tests for core functionality.
"""

import pytest
import numpy as np
import pandas as pd
from scmetabolism import ScMetabolism
from scmetabolism.utils import load_gene_sets, alra_imputation


class TestScMetabolism:
    """Test ScMetabolism class."""
    
    def setup_method(self):
        """Set up test data."""
        # Create synthetic data
        np.random.seed(42)
        self.n_genes = 100
        self.n_cells = 50
        
        self.genes = [f'GENE_{i}' for i in range(self.n_genes)]
        self.cells = [f'CELL_{i}' for i in range(self.n_cells)]
        
        # Create count matrix
        self.count_matrix = pd.DataFrame(
            np.random.negative_binomial(5, 0.3, size=(self.n_genes, self.n_cells)),
            index=self.genes,
            columns=self.cells
        )
        
        # Create ScMetabolism instance
        self.sc_metab = ScMetabolism()
    
    def test_initialization(self):
        """Test ScMetabolism initialization."""
        assert self.sc_metab.metabolism_scores is None
        assert self.sc_metab.gene_sets is None
        assert self.sc_metab.method_used is None
    
    def test_compute_metabolism_aucell(self):
        """Test metabolism computation with AUCell."""
        # Mock gene sets for testing
        mock_gene_sets = {
            'Pathway1': ['GENE_1', 'GENE_2', 'GENE_3'],
            'Pathway2': ['GENE_4', 'GENE_5', 'GENE_6'],
            'Pathway3': ['GENE_7', 'GENE_8', 'GENE_9']
        }
        
        # Set mock gene sets
        self.sc_metab.gene_sets = mock_gene_sets
        
        # Test AUCell computation
        scores = self.sc_metab._compute_aucell(self.count_matrix, n_cores=1)
        
        assert isinstance(scores, pd.DataFrame)
        assert scores.shape[0] == len(mock_gene_sets)
        assert scores.shape[1] == self.n_cells
        assert all(pathway in scores.index for pathway in mock_gene_sets.keys())
    
    def test_compute_metabolism_ssgsea(self):
        """Test metabolism computation with ssGSEA."""
        # Mock gene sets
        mock_gene_sets = {
            'Pathway1': ['GENE_1', 'GENE_2', 'GENE_3'],
            'Pathway2': ['GENE_4', 'GENE_5', 'GENE_6']
        }
        
        self.sc_metab.gene_sets = mock_gene_sets
        
        # Test ssGSEA computation
        scores = self.sc_metab._compute_ssgsea(self.count_matrix, n_cores=1)
        
        assert isinstance(scores, pd.DataFrame)
        assert scores.shape[0] == len(mock_gene_sets)
        assert scores.shape[1] == self.n_cells
    
    def test_compute_metabolism_gsva(self):
        """Test metabolism computation with GSVA."""
        # Mock gene sets
        mock_gene_sets = {
            'Pathway1': ['GENE_1', 'GENE_2', 'GENE_3'],
            'Pathway2': ['GENE_4', 'GENE_5', 'GENE_6']
        }
        
        self.sc_metab.gene_sets = mock_gene_sets
        
        # Test GSVA computation
        scores = self.sc_metab._compute_gsva(self.count_matrix, n_cores=1)
        
        assert isinstance(scores, pd.DataFrame)
        assert scores.shape[0] == len(mock_gene_sets)
        assert scores.shape[1] == self.n_cells
    
    def test_invalid_method(self):
        """Test invalid method raises error."""
        with pytest.raises(ValueError, match="Unknown method"):
            self.sc_metab.compute_metabolism(
                self.count_matrix,
                method="invalid_method"
            )


class TestUtils:
    """Test utility functions."""
    
    def test_alra_imputation(self):
        """Test ALRA imputation."""
        # Create test data with some zeros
        np.random.seed(42)
        data = pd.DataFrame(
            np.random.negative_binomial(2, 0.5, size=(50, 30)),
            index=[f'Gene_{i}' for i in range(50)],
            columns=[f'Cell_{i}' for i in range(30)]
        )
        
        # Add some zeros
        data.iloc[:10, :5] = 0
        
        # Test imputation
        imputed = alra_imputation(data, k=10)
        
        assert isinstance(imputed, pd.DataFrame)
        assert imputed.shape == data.shape
        assert imputed.index.equals(data.index)
        assert imputed.columns.equals(data.columns)
        
        # Check that zeros are preserved where original was zero
        original_zeros = (data == 0)
        assert (imputed[original_zeros] == 0).all().all()


if __name__ == "__main__":
    pytest.main([__file__])