"""
Tests for QC Metrics Calculation

These tests validate that QC metrics (library size, detected genes,
mitochondrial content) are calculated correctly.
"""

import pytest
import pandas as pd
import numpy as np


class TestLibrarySize:
    """Test library size calculation."""
    
    def test_library_size_sum(self, dataset_from_simple):
        """
        Library size should be the sum of all gene counts per sample.
        
        Sample1: 100 + 300 + 600 + 0 = 1000
        Sample2: 200 + 400 + 400 + 0 = 1000
        Sample3: 50 + 100 + 150 + 700 = 1000
        """
        dataset = dataset_from_simple
        dataset.calculate_qc_metrics()
        
        assert dataset.qc_metrics.loc['Sample1', 'Library_Size'] == 1000
        assert dataset.qc_metrics.loc['Sample2', 'Library_Size'] == 1000
        assert dataset.qc_metrics.loc['Sample3', 'Library_Size'] == 1000
    
    def test_library_size_unequal(self, dataset_from_uneven):
        """
        Library size should correctly calculate different values per sample.
        
        Sample1: 100 + 100 = 200
        Sample2: 1000 + 1000 = 2000
        Sample3: 10 + 10 = 20
        """
        dataset = dataset_from_uneven
        dataset.calculate_qc_metrics()
        
        assert dataset.qc_metrics.loc['Sample1', 'Library_Size'] == 200
        assert dataset.qc_metrics.loc['Sample2', 'Library_Size'] == 2000
        assert dataset.qc_metrics.loc['Sample3', 'Library_Size'] == 20


class TestDetectedGenes:
    """Test detected genes calculation."""
    
    def test_detected_genes_count(self, dataset_from_simple):
        """
        Detected genes should count genes with >0 expression.
        
        Sample1: GeneA(100), GeneB(300), GeneC(600), MT-ATP6(0) → 3 detected
        Sample2: GeneA(200), GeneB(400), GeneC(400), MT-ATP6(0) → 3 detected
        Sample3: GeneA(50), GeneB(100), GeneC(150), MT-ATP6(700) → 4 detected
        """
        dataset = dataset_from_simple
        dataset.calculate_qc_metrics()
        
        assert dataset.qc_metrics.loc['Sample1', 'Detected_Genes'] == 3
        assert dataset.qc_metrics.loc['Sample2', 'Detected_Genes'] == 3
        assert dataset.qc_metrics.loc['Sample3', 'Detected_Genes'] == 4
    
    def test_detected_genes_sparse(self, dataset_from_sparse):
        """
        Test detected genes with sparse matrix.
        
        Sample1: GeneA(10), GeneB(0), GeneC(1), GeneD(0) → 2 detected
        Sample2: GeneA(20), GeneB(0), GeneC(0), GeneD(0) → 1 detected
        Sample3: GeneA(30), GeneB(50), GeneC(2), GeneD(0) → 3 detected
        Sample4: GeneA(40), GeneB(60), GeneC(0), GeneD(0) → 2 detected
        """
        dataset = dataset_from_sparse
        dataset.calculate_qc_metrics()
        
        assert dataset.qc_metrics.loc['Sample1', 'Detected_Genes'] == 2
        assert dataset.qc_metrics.loc['Sample2', 'Detected_Genes'] == 1
        assert dataset.qc_metrics.loc['Sample3', 'Detected_Genes'] == 3
        assert dataset.qc_metrics.loc['Sample4', 'Detected_Genes'] == 2


class TestMitochondrialContent:
    """Test mitochondrial content calculation."""
    
    def test_mt_percentage_calculation(self, dataset_from_simple):
        """
        MT% should be the percentage of reads from MT-* genes.
        
        Sample1: MT counts = 0, total = 1000 → MT% = 0%
        Sample2: MT counts = 0, total = 1000 → MT% = 0%
        Sample3: MT counts = 700, total = 1000 → MT% = 70%
        """
        dataset = dataset_from_simple
        dataset.calculate_qc_metrics()
        
        assert dataset.qc_metrics.loc['Sample1', 'MT_Percentage'] == pytest.approx(0.0)
        assert dataset.qc_metrics.loc['Sample2', 'MT_Percentage'] == pytest.approx(0.0)
        assert dataset.qc_metrics.loc['Sample3', 'MT_Percentage'] == pytest.approx(70.0)
    
    def test_mt_no_mt_genes_returns_zero(self, dataset_from_uneven):
        """
        If no MT-* genes exist, MT% should be 0 for all samples.
        """
        dataset = dataset_from_uneven
        # This matrix has no MT-* genes
        dataset.calculate_qc_metrics()
        
        for sample in dataset.sample_ids:
            assert dataset.qc_metrics.loc[sample, 'MT_Percentage'] == 0.0


class TestOutlierDetection:
    """Test outlier detection functionality."""
    
    def test_flag_low_library_size(self, dataset_from_uneven):
        """
        Samples below library size threshold should be flagged.
        
        Sample1: 200 reads
        Sample2: 2000 reads
        Sample3: 20 reads
        
        With threshold = 100:
        - Sample3 should fail (20 < 100)
        """
        dataset = dataset_from_uneven
        dataset.calculate_qc_metrics()
        
        failed = dataset.flag_outliers(lib_size_min=100)
        
        assert 'Sample3' in failed.index
        assert 'Sample1' not in failed.index  # 200 >= 100
        assert 'Sample2' not in failed.index  # 2000 >= 100
    
    def test_flag_high_mt_percentage(self, dataset_from_simple):
        """
        Samples above MT% threshold should be flagged.
        
        Sample1: MT% = 0%
        Sample2: MT% = 0%
        Sample3: MT% = 70%
        
        With threshold = 20%:
        - Sample3 should fail (70% > 20%)
        """
        dataset = dataset_from_simple
        dataset.calculate_qc_metrics()
        
        failed = dataset.flag_outliers(lib_size_min=0, mt_threshold=20.0)
        
        assert 'Sample3' in failed.index
        assert 'Sample1' not in failed.index
        assert 'Sample2' not in failed.index
    
    def test_multiple_failure_reasons(self, dataset_from_simple):
        """
        Samples can fail for multiple reasons.
        """
        # Create a dataset where a sample fails both criteria
        df = pd.DataFrame({
            'GeneA': [100, 10],
            'MT-X': [0, 90]
        }, index=['Good', 'Bad'])
        
        from src.dataset import RNASeqDataset
        dataset = RNASeqDataset.from_dataframe(df)
        dataset.calculate_qc_metrics()
        
        # 'Bad' has library size 100, MT% = 90%
        failed = dataset.flag_outliers(lib_size_min=500, mt_threshold=20.0)
        
        assert 'Bad' in failed.index
        # Check both reasons are listed
        reasons = failed.loc['Bad', 'Failure_Reasons']
        assert 'Low library size' in reasons
        assert 'High MT%' in reasons
