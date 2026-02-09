"""
Pytest Configuration and Fixtures

This module provides shared test fixtures for RNA-Seq pipeline testing,
including pre-defined small matrices where results can be calculated manually.
"""

import pytest
import pandas as pd
import numpy as np
import tempfile
import os
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.dataset import RNASeqDataset


@pytest.fixture
def simple_expression_matrix():
    """
    A simple 3-sample × 4-gene expression matrix for manual verification.
    
    Matrix:
                GeneA  GeneB  GeneC  MT-ATP6
    Sample1      100    300    600      0
    Sample2      200    400    400      0
    Sample3       50    100    150    700
    
    Expected totals:
    - Sample1: 1000 reads, MT% = 0%
    - Sample2: 1000 reads, MT% = 0%
    - Sample3: 1000 reads, MT% = 70%
    
    Expected CPM (all samples have same library size = 1000):
    - CPM = count * 1,000,000 / 1000 = count * 1000
    """
    return pd.DataFrame({
        'GeneA': [100, 200, 50],
        'GeneB': [300, 400, 100],
        'GeneC': [600, 400, 150],
        'MT-ATP6': [0, 0, 700]  # Mitochondrial gene
    }, index=['Sample1', 'Sample2', 'Sample3'])


@pytest.fixture
def uneven_library_matrix():
    """
    Expression matrix with different library sizes to test normalization.
    
    Matrix:
                GeneA  GeneB
    Sample1      100    100   (total = 200)
    Sample2     1000   1000   (total = 2000)
    Sample3       10     10   (total = 20)
    
    Expected CPM:
    - Sample1: [500000, 500000]
    - Sample2: [500000, 500000]
    - Sample3: [500000, 500000]
    
    After normalization, all samples should have the same CPM values
    because they have the same proportions.
    """
    return pd.DataFrame({
        'GeneA': [100, 1000, 10],
        'GeneB': [100, 1000, 10]
    }, index=['Sample1', 'Sample2', 'Sample3'])


@pytest.fixture
def sparse_expression_matrix():
    """
    Expression matrix with many zero values for testing gene filtering.
    
    Matrix:
                GeneA  GeneB  GeneC  GeneD
    Sample1       10      0      1      0
    Sample2       20      0      0      0
    Sample3       30     50      2      0
    Sample4       40     60      0      0
    
    With min_cpm=1.0 and min_samples_pct=0.5 (i.e., 2 samples):
    - GeneA: expressed in 4/4 samples ✓
    - GeneB: expressed in 2/4 samples ✓
    - GeneC: expressed in 2/4 samples ✓
    - GeneD: expressed in 0/4 samples ✗
    """
    return pd.DataFrame({
        'GeneA': [10, 20, 30, 40],
        'GeneB': [0, 0, 50, 60],
        'GeneC': [1, 0, 2, 0],
        'GeneD': [0, 0, 0, 0]
    }, index=['Sample1', 'Sample2', 'Sample3', 'Sample4'])


@pytest.fixture
def dataset_from_simple(simple_expression_matrix):
    """Create an RNASeqDataset from the simple matrix."""
    return RNASeqDataset.from_dataframe(simple_expression_matrix)


@pytest.fixture
def dataset_from_uneven(uneven_library_matrix):
    """Create an RNASeqDataset from the uneven library matrix."""
    return RNASeqDataset.from_dataframe(uneven_library_matrix)


@pytest.fixture
def dataset_from_sparse(sparse_expression_matrix):
    """Create an RNASeqDataset from the sparse matrix."""
    return RNASeqDataset.from_dataframe(sparse_expression_matrix)


@pytest.fixture
def temp_csv_file(simple_expression_matrix, tmp_path):
    """
    Create a temporary CSV file in the expected input format (genes × samples).
    This mimics the TCGA format where genes are rows and samples are columns.
    """
    # Transpose to get genes × samples (original file format)
    genes_x_samples = simple_expression_matrix.T
    genes_x_samples.index.name = "GeneID"
    
    filepath = tmp_path / "test_expression.csv"
    genes_x_samples.to_csv(filepath)
    
    return str(filepath)


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary output directory."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return str(output_dir)
