"""
Tests for CPM Normalization

These tests validate that the CPM (Counts Per Million) normalization
produces mathematically correct results.

CPM Formula: CPM = (count / total_counts) × 10^6
"""

import pytest
import pandas as pd
import numpy as np


class TestCPMNormalization:
    """Test suite for CPM normalization functionality."""
    
    def test_cpm_equal_library_sizes(self, dataset_from_simple):
        """
        When all samples have equal library sizes, CPM should scale
        counts by a constant factor.
        
        Given: All samples have 1000 total reads
        Expected: CPM = count * (1,000,000 / 1,000) = count * 1000
        """
        dataset = dataset_from_simple
        dataset.normalize(method="cpm")
        
        # Sample1, GeneA: raw = 100, expected CPM = 100 * 1000 = 100,000
        assert dataset.normalized.loc['Sample1', 'GeneA'] == pytest.approx(100000)
        
        # Sample2, GeneB: raw = 400, expected CPM = 400 * 1000 = 400,000
        assert dataset.normalized.loc['Sample2', 'GeneB'] == pytest.approx(400000)
        
        # Sample3, MT-ATP6: raw = 700, expected CPM = 700 * 1000 = 700,000
        assert dataset.normalized.loc['Sample3', 'MT-ATP6'] == pytest.approx(700000)
    
    def test_cpm_unequal_library_sizes(self, dataset_from_uneven):
        """
        When samples have different library sizes, CPM should normalize
        them to be comparable.
        
        Given: Samples have 200, 2000, 20 total reads respectively
               but all have 50/50 split between GeneA and GeneB
        Expected: All samples should have CPM of 500,000 for each gene
        """
        dataset = dataset_from_uneven
        dataset.normalize(method="cpm")
        
        # All samples should have 50% GeneA and 50% GeneB
        # So CPM should be 500,000 for each gene in each sample
        for sample in ['Sample1', 'Sample2', 'Sample3']:
            assert dataset.normalized.loc[sample, 'GeneA'] == pytest.approx(500000)
            assert dataset.normalized.loc[sample, 'GeneB'] == pytest.approx(500000)
    
    def test_cpm_row_sums_equal_million(self, dataset_from_simple):
        """
        After CPM normalization, each sample (row) should sum to 1,000,000.
        """
        dataset = dataset_from_simple
        dataset.normalize(method="cpm")
        
        row_sums = dataset.normalized.sum(axis=1)
        
        for sample in dataset.sample_ids:
            assert row_sums[sample] == pytest.approx(1_000_000)
    
    def test_cpm_preserves_proportions(self, dataset_from_simple):
        """
        CPM normalization should preserve the relative proportions
        of genes within each sample.
        """
        dataset = dataset_from_simple
        
        # Calculate raw proportions for Sample1
        raw_proportions = dataset.raw_counts.loc['Sample1'] / dataset.raw_counts.loc['Sample1'].sum()
        
        # Normalize
        dataset.normalize(method="cpm")
        
        # Calculate normalized proportions (divide by 1M since CPM)
        norm_proportions = dataset.normalized.loc['Sample1'] / 1_000_000
        
        # Raw and normalized proportions should be equal
        pd.testing.assert_series_equal(
            raw_proportions, 
            norm_proportions,
            check_names=False
        )
    
    def test_normalize_requires_loaded_data(self):
        """Normalizing without loading data should raise an error."""
        from src.dataset import RNASeqDataset
        
        dataset = RNASeqDataset()  # No data loaded
        
        with pytest.raises(ValueError, match="No data loaded"):
            dataset.normalize()
    
    def test_invalid_normalization_method(self, dataset_from_simple):
        """An invalid normalization method should raise an error."""
        dataset = dataset_from_simple
        
        with pytest.raises(ValueError, match="Unknown normalization method"):
            dataset.normalize(method="invalid_method")


class TestNormalizationEdgeCases:
    """Test edge cases for normalization."""
    
    def test_zero_counts_remain_zero(self, dataset_from_simple):
        """Genes with zero counts should remain zero after normalization."""
        dataset = dataset_from_simple
        dataset.normalize(method="cpm")
        
        # Sample1 and Sample2 have 0 counts for MT-ATP6
        assert dataset.normalized.loc['Sample1', 'MT-ATP6'] == 0
        assert dataset.normalized.loc['Sample2', 'MT-ATP6'] == 0
    
    def test_normalization_creates_new_dataframe(self, dataset_from_simple):
        """
        Normalization should not modify the raw counts.
        """
        dataset = dataset_from_simple
        raw_before = dataset.raw_counts.copy()
        
        dataset.normalize(method="cpm")
        
        pd.testing.assert_frame_equal(dataset.raw_counts, raw_before)
