"""
Performance optimization utilities for scMetabolism.
"""

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse
from numba import jit, prange
import warnings
from typing import Union, Optional
import time
from functools import wraps


def timer(func):
    """Decorator to time function execution."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} completed in {end_time - start_time:.2f} seconds")
        return result
    return wrapper


@jit(nopython=True, parallel=True)
def fast_ranking(matrix):
    """Fast ranking using numba for large matrices."""
    n_genes, n_cells = matrix.shape
    rankings = np.zeros((n_genes, n_cells), dtype=np.float32)
    
    for cell_idx in prange(n_cells):
        cell_data = matrix[:, cell_idx]
        # Get sorted indices
        sorted_indices = np.argsort(-cell_data)  # Descending order
        
        # Assign ranks
        for rank, gene_idx in enumerate(sorted_indices):
            rankings[gene_idx, cell_idx] = rank + 1
            
    return rankings


class SparseMatrixHandler:
    """Handle sparse matrices efficiently."""
    
    @staticmethod
    def to_sparse(matrix: Union[pd.DataFrame, np.ndarray], 
                  threshold: float = 0.7) -> Union[csr_matrix, pd.DataFrame]:
        """
        Convert to sparse matrix if sparsity is above threshold.
        
        Parameters:
        -----------
        matrix : pd.DataFrame or np.ndarray
            Input matrix
        threshold : float
            Sparsity threshold (fraction of zeros)
            
        Returns:
        --------
        Sparse or dense matrix
        """
        if isinstance(matrix, pd.DataFrame):
            values = matrix.values
        else:
            values = matrix
            
        sparsity = np.count_nonzero(values == 0) / values.size
        
        if sparsity > threshold:
            print(f"Converting to sparse matrix (sparsity: {sparsity:.2f})")
            if isinstance(matrix, pd.DataFrame):
                return csr_matrix(values), matrix.index, matrix.columns
            else:
                return csr_matrix(values)
        
        return matrix
    
    @staticmethod
    def sparse_ranking(sparse_matrix: csr_matrix) -> np.ndarray:
        """Efficient ranking for sparse matrices."""
        dense_matrix = sparse_matrix.toarray()
        return fast_ranking(dense_matrix.astype(np.float32))


class MemoryOptimizer:
    """Memory optimization utilities."""
    
    @staticmethod
    def optimize_dtypes(df: pd.DataFrame) -> pd.DataFrame:
        """Optimize DataFrame dtypes to reduce memory usage."""
        original_memory = df.memory_usage(deep=True).sum()
        
        # Optimize numeric columns
        for col in df.select_dtypes(include=[np.number]).columns:
            col_min = df[col].min()
            col_max = df[col].max()
            
            if df[col].dtype == np.int64:
                if col_min >= np.iinfo(np.int8).min and col_max <= np.iinfo(np.int8).max:
                    df[col] = df[col].astype(np.int8)
                elif col_min >= np.iinfo(np.int16).min and col_max <= np.iinfo(np.int16).max:
                    df[col] = df[col].astype(np.int16)
                elif col_min >= np.iinfo(np.int32).min and col_max <= np.iinfo(np.int32).max:
                    df[col] = df[col].astype(np.int32)
            
            elif df[col].dtype == np.float64:
                if col_min >= np.finfo(np.float32).min and col_max <= np.finfo(np.float32).max:
                    df[col] = df[col].astype(np.float32)
        
        # Optimize object columns
        for col in df.select_dtypes(include=['object']).columns:
            if df[col].nunique() / len(df) < 0.5:  # If less than 50% unique values
                df[col] = df[col].astype('category')
        
        optimized_memory = df.memory_usage(deep=True).sum()
        reduction = (original_memory - optimized_memory) / original_memory * 100
        
        print(f"Memory usage reduced by {reduction:.1f}% "
              f"({original_memory/1024**2:.1f}MB -> {optimized_memory/1024**2:.1f}MB)")
        
        return df
    
    @staticmethod
    def chunk_processor(data: pd.DataFrame, chunk_size: int = 1000):
        """Process data in chunks to reduce memory usage."""
        n_chunks = len(data) // chunk_size + (1 if len(data) % chunk_size else 0)
        
        for i in range(n_chunks):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, len(data))
            yield data.iloc[start_idx:end_idx]


class BatchProcessor:
    """Process large datasets in batches."""
    
    def __init__(self, batch_size: int = 5000):
        self.batch_size = batch_size
    
    def process_cells_in_batches(self, count_matrix: pd.DataFrame, 
                                processing_func, **kwargs):
        """Process cells in batches."""
        n_cells = count_matrix.shape[1]
        n_batches = n_cells // self.batch_size + (1 if n_cells % self.batch_size else 0)
        
        results = []
        
        for batch_idx in range(n_batches):
            start_idx = batch_idx * self.batch_size
            end_idx = min((batch_idx + 1) * self.batch_size, n_cells)
            
            batch_data = count_matrix.iloc[:, start_idx:end_idx]
            batch_result = processing_func(batch_data, **kwargs)
            results.append(batch_result)
            
            print(f"Processed batch {batch_idx + 1}/{n_batches}")
        
        # Combine results
        if isinstance(results[0], pd.DataFrame):
            return pd.concat(results, axis=1)
        else:
            return np.concatenate(results, axis=1)


def profile_memory_usage(func):
    """Decorator to profile memory usage."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            import psutil
            import os
            
            process = psutil.Process(os.getpid())
            mem_before = process.memory_info().rss / 1024 / 1024  # MB
            
            result = func(*args, **kwargs)
            
            mem_after = process.memory_info().rss / 1024 / 1024  # MB
            mem_diff = mem_after - mem_before
            
            print(f"{func.__name__} memory usage: {mem_diff:.1f}MB "
                  f"(before: {mem_before:.1f}MB, after: {mem_after:.1f}MB)")
            
            return result
            
        except ImportError:
            warnings.warn("psutil not available for memory profiling")
            return func(*args, **kwargs)
    
    return wrapper