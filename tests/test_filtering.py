"""
Tests for Gene Filtering

These tests validate that the gene filtering logic correctly removes
lowly-expressed genes based on CPM thresholds.

Filter Rule: Keep genes with >min_cpm in at least min_samples_pct of samples
"""

import pytest
import pandas as pd
import numpy as np


class TestGeneFiltering:
    """Test suite for gene filtering functionality."""
    
    def test_filter_removes_unexpressed_genes(self, dataset_from_sparse):
        """
        Genes with zero expression across all samples should be filtered out.
        
        GeneD has 0 counts in all 4 samples → should be removed
        """
        dataset = dataset_from_sparse
        dataset.normalize(method="cpm")
        filtered = dataset.filter_genes(min_cpm=1.0, min_samples_pct=0.5)
        
        assert 'GeneD' not in filtered
    
    def test_filter_keeps_high_expression_genes(self, dataset_from_sparse):
        """
        Genes with high expression in most samples should be kept.
        
        GeneA is expressed in 4/4 samples → should be kept
        """
        dataset = dataset_from_sparse
        dataset.normalize(method="cpm")
        filtered = dataset.filter_genes(min_cpm=1.0, min_samples_pct=0.5)
        
        assert 'GeneA' in filtered
    
    def test_filter_respects_sample_percentage(self, dataset_from_sparse):
        """
        Genes must be expressed in at least min_samples_pct of samples.
        
        With min_samples_pct=0.5 (2/4 samples):
        - GeneB: expressed in 2/4 samples ✓
        
        With min_samples_pct=0.75 (3/4 samples):
        - GeneB: expressed in 2/4 samples ✗
        """
        dataset = dataset_from_sparse
        dataset.normalize(method="cpm")
        
        # 50% threshold - GeneB should pass
        filtered_50 = dataset.filter_genes(min_cpm=1.0, min_samples_pct=0.5)
        assert 'GeneB' in filtered_50
        
        # Reset and apply 75% threshold
        dataset.normalized = None
        dataset.normalize(method="cpm")
        filtered_75 = dataset.filter_genes(min_cpm=1.0, min_samples_pct=0.75)
        assert 'GeneB' not in filtered_75
    
    def test_filter_respects_cpm_threshold(self, dataset_from_sparse):
        """
        Only genes exceeding the min_cpm threshold should count as expressed.
        """
        dataset = dataset_from_sparse
        dataset.normalize(method="cpm")
        
        # Very high CPM threshold should filter out more genes
        filtered_high = dataset.filter_genes(min_cpm=100000, min_samples_pct=0.5)
        
        # With such a high threshold, most genes won't pass
        assert len(filtered_high) < 4
    
    def test_get_filtered_data_returns_subset(self, dataset_from_sparse):
        """
        get_filtered_data() should return data for filtered genes only.
        """
        dataset = dataset_from_sparse
        dataset.normalize(method="cpm")
        dataset.filter_genes(min_cpm=1.0, min_samples_pct=0.5)
        
        filtered_data = dataset.get_filtered_data(normalized=True)
        
        # Should only contain columns for filtered genes
        assert set(filtered_data.columns) == set(dataset.filtered_genes)
        
        # Should still have all samples (rows)
        assert len(filtered_data) == len(dataset.sample_ids)
    
    def test_filter_auto_normalizes_if_needed(self, dataset_from_sparse):
        """
        If normalization hasn't been run, filter_genes should run it automatically.
        """
        dataset = dataset_from_sparse
        # Don't normalize explicitly
        
        # This should work and auto-normalize
        filtered = dataset.filter_genes(min_cpm=1.0, min_samples_pct=0.5)
        
        assert dataset.normalized is not None
        assert len(filtered) > 0


class TestFilteringEdgeCases:
    """Test edge cases for filtering."""
    
    def test_filter_100_percent_threshold(self, dataset_from_simple):
        """
        With 100% sample requirement, only genes expressed in ALL samples pass.
        """
        dataset = dataset_from_simple
        dataset.normalize(method="cpm")
        
        # MT-ATP6 is only expressed in Sample3
        filtered = dataset.filter_genes(min_cpm=1.0, min_samples_pct=1.0)
        
        assert 'MT-ATP6' not in filtered
        # GeneA, GeneB, GeneC are expressed in all 3 samples
        assert 'GeneA' in filtered
        assert 'GeneB' in filtered
        assert 'GeneC' in filtered
    
    def test_filter_zero_percent_threshold(self, dataset_from_simple):
        """
        With 0% sample requirement, all genes should pass.
        """
        dataset = dataset_from_simple
        dataset.normalize(method="cpm")
        
        filtered = dataset.filter_genes(min_cpm=1.0, min_samples_pct=0.0)
        
        # All 4 genes should pass
        assert len(filtered) == 4
