"""
Tests for Data Loading

These tests validate that data loading from CSV files works correctly.
"""

import pytest
import pandas as pd
import numpy as np
import os


class TestDataLoading:
    """Test data loading functionality."""
    
    def test_load_from_csv(self, temp_csv_file):
        """
        Loading from CSV should correctly transpose and set up the matrix.
        """
        from src.dataset import RNASeqDataset
        
        dataset = RNASeqDataset(filepath=temp_csv_file)
        
        # Should have 3 samples (from simple_expression_matrix)
        assert len(dataset.sample_ids) == 3
        
        # Should have 4 genes
        assert len(dataset.gene_ids) == 4
        
        # Data should be transposed (samples × genes)
        assert dataset.raw_counts.shape == (3, 4)
    
    def test_load_preserves_values(self, temp_csv_file, simple_expression_matrix):
        """
        Loaded values should match the original matrix.
        """
        from src.dataset import RNASeqDataset
        
        dataset = RNASeqDataset(filepath=temp_csv_file)
        
        # Check specific values
        assert dataset.raw_counts.loc['Sample1', 'GeneA'] == 100
        assert dataset.raw_counts.loc['Sample3', 'MT-ATP6'] == 700
    
    def test_load_nonexistent_file_raises_error(self):
        """
        Loading a nonexistent file should raise FileNotFoundError.
        """
        from src.dataset import RNASeqDataset
        
        with pytest.raises(FileNotFoundError):
            RNASeqDataset(filepath="nonexistent_file.csv")
    
    def test_from_dataframe(self, simple_expression_matrix):
        """
        from_dataframe should create a dataset from an existing DataFrame.
        """
        from src.dataset import RNASeqDataset
        
        dataset = RNASeqDataset.from_dataframe(simple_expression_matrix)
        
        assert len(dataset.sample_ids) == 3
        assert len(dataset.gene_ids) == 4
        pd.testing.assert_frame_equal(dataset.raw_counts, simple_expression_matrix)


class TestDatasetRepresentation:
    """Test dataset string representation and summary."""
    
    def test_repr(self, dataset_from_simple):
        """
        __repr__ should return a descriptive string.
        """
        repr_str = repr(dataset_from_simple)
        
        assert "RNASeqDataset" in repr_str
        assert "3 samples" in repr_str
        assert "4 genes" in repr_str
    
    def test_get_summary(self, dataset_from_simple):
        """
        get_summary() should return a dictionary with dataset info.
        """
        dataset = dataset_from_simple
        dataset.normalize()
        dataset.calculate_qc_metrics()
        
        summary = dataset.get_summary()
        
        assert summary["n_samples"] == 3
        assert summary["n_genes_total"] == 4
        assert summary["normalized"] == True
        assert summary["qc_calculated"] == True
