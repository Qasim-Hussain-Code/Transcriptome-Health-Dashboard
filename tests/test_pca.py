"""
Tests for PCA Analysis

These tests validate that PCA dimensionality reduction works correctly.
"""

import pytest
import pandas as pd
import numpy as np


class TestPCA:
    """Test PCA functionality."""
    
    def test_pca_returns_correct_shape(self, dataset_from_simple):
        """
        PCA should return coordinates with shape (n_samples, n_components).
        """
        dataset = dataset_from_simple
        dataset.normalize(method="cpm")
        
        coords, variance = dataset.run_pca(n_components=2)
        
        # Should have 3 samples and 2 components
        assert coords.shape == (3, 2)
        
        # Should have 2 variance values
        assert len(variance) == 2
    
    def test_pca_variance_sums_to_less_than_one(self, dataset_from_simple):
        """
        Explained variance ratios should sum to ≤ 1.0.
        """
        dataset = dataset_from_simple
        dataset.normalize(method="cpm")
        
        _, variance = dataset.run_pca(n_components=2)
        
        assert sum(variance) <= 1.0
        assert all(v >= 0 for v in variance)  # All non-negative
    
    def test_pca_stores_results(self, dataset_from_simple):
        """
        PCA results should be stored in dataset.pca_results.
        """
        dataset = dataset_from_simple
        dataset.normalize(method="cpm")
        
        dataset.run_pca(n_components=2)
        
        assert dataset.pca_results is not None
        assert "coordinates" in dataset.pca_results
        assert "explained_variance" in dataset.pca_results
        assert "n_components" in dataset.pca_results
        assert dataset.pca_results["n_components"] == 2
    
    def test_pca_uses_filtered_genes(self, dataset_from_sparse):
        """
        When use_filtered=True, PCA should only use filtered genes.
        """
        dataset = dataset_from_sparse
        dataset.normalize(method="cpm")
        dataset.filter_genes(min_cpm=1.0, min_samples_pct=0.5)
        
        coords, _ = dataset.run_pca(n_components=2, use_filtered=True)
        
        # Should still have all samples
        assert coords.shape[0] == len(dataset.sample_ids)
    
    def test_pca_auto_normalizes(self, dataset_from_simple):
        """
        PCA should auto-normalize if not already done.
        """
        dataset = dataset_from_simple
        # Don't normalize explicitly
        
        coords, _ = dataset.run_pca(n_components=2, use_filtered=False)
        
        assert dataset.normalized is not None
        assert coords.shape == (3, 2)


class TestPCAEdgeCases:
    """Test PCA edge cases."""
    
    def test_pca_single_component(self, dataset_from_simple):
        """
        PCA should work with n_components=1.
        """
        dataset = dataset_from_simple
        dataset.normalize(method="cpm")
        
        coords, variance = dataset.run_pca(n_components=1)
        
        assert coords.shape == (3, 1)
        assert len(variance) == 1
    
    def test_pca_max_components(self, dataset_from_simple):
        """
        Number of components is limited by min(n_samples, n_features).
        """
        dataset = dataset_from_simple  # 3 samples × 4 genes
        dataset.normalize(method="cpm")
        
        # Request 3 components (max possible for 3 samples)
        coords, variance = dataset.run_pca(n_components=3)
        
        assert coords.shape == (3, 3)
