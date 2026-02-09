"""
Core Dataset Class for RNA-Seq Analysis

This module provides the RNASeqDataset class that encapsulates all data operations
for RNA-Seq quality control and analysis, including:
- Data loading and preprocessing
- CPM normalization
- Gene filtering
- QC metrics calculation (library size, detected genes, MT content)
- PCA dimensionality reduction
- Outlier detection and flagging
"""

import pandas as pd
import numpy as np
import logging
import time
import os
from typing import Optional, Tuple, List, Dict, Any
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


class RNASeqDataset:
    """
    A comprehensive container for RNA-Seq expression data with built-in
    quality control, normalization, and analysis methods.
    
    Attributes:
        raw_counts (pd.DataFrame): Original counts matrix (samples × genes)
        normalized (pd.DataFrame): CPM-normalized expression data
        qc_metrics (pd.DataFrame): Per-sample quality control statistics
        filtered_genes (List[str]): List of genes passing expression filter
        sample_ids (List[str]): List of sample/patient identifiers
        gene_ids (List[str]): List of gene identifiers
        logger (logging.Logger): Logger instance for pipeline tracking
    
    Example:
        >>> dataset = RNASeqDataset("data/expression.csv")
        >>> dataset.normalize(method="cpm")
        >>> dataset.filter_genes(min_cpm=1.0, min_samples_pct=0.5)
        >>> dataset.calculate_qc_metrics()
        >>> outliers = dataset.flag_outliers(lib_threshold=20_000_000)
    """
    
    def __init__(
        self,
        filepath: Optional[str] = None,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the RNASeqDataset.
        
        Args:
            filepath: Path to CSV file containing expression matrix
                      (genes × samples format, first column = GeneID)
            logger: Optional logger instance. If None, creates a default logger.
        """
        # Initialize logger
        self.logger = logger or self._create_default_logger()
        
        # Data containers
        self.raw_counts: Optional[pd.DataFrame] = None
        self.normalized: Optional[pd.DataFrame] = None
        self.qc_metrics: Optional[pd.DataFrame] = None
        self.pca_results: Optional[Dict[str, Any]] = None
        
        # Metadata
        self.filtered_genes: List[str] = []
        self.sample_ids: List[str] = []
        self.gene_ids: List[str] = []
        self.filepath: Optional[str] = filepath
        
        # Load data if filepath provided
        if filepath:
            self.load_data(filepath)
    
    @classmethod
    def from_dataframe(
        cls,
        df: pd.DataFrame,
        logger: Optional[logging.Logger] = None
    ) -> 'RNASeqDataset':
        """
        Create RNASeqDataset from an existing DataFrame.
        
        Args:
            df: Expression DataFrame (samples × genes)
            logger: Optional logger instance
            
        Returns:
            Initialized RNASeqDataset instance
        """
        instance = cls(filepath=None, logger=logger)
        instance.raw_counts = df.copy()
        instance.sample_ids = list(df.index)
        instance.gene_ids = list(df.columns)
        instance.logger.info(
            f"Created dataset from DataFrame: {len(instance.sample_ids)} samples × "
            f"{len(instance.gene_ids)} genes"
        )
        return instance
    
    def _create_default_logger(self) -> logging.Logger:
        """Create a default logger with console output."""
        logger = logging.getLogger("RNASeqDataset")
        if not logger.handlers:
            logger.setLevel(logging.INFO)
            handler = logging.StreamHandler()
            handler.setFormatter(
                logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            )
            logger.addHandler(handler)
        return logger
    
    def load_data(self, filepath: str) -> pd.DataFrame:
        """
        Load and preprocess expression data from CSV file.
        
        The expected format is genes × samples with first column as GeneID.
        Data is automatically transposed to samples × genes format.
        
        Args:
            filepath: Path to the CSV file
            
        Returns:
            Transposed DataFrame (samples × genes)
            
        Raises:
            FileNotFoundError: If the file doesn't exist
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Error: The file '{filepath}' was not found.")
        
        self.logger.info(f"Loading expression data: {filepath}")
        start_time = time.time()
        
        # Load raw CSV
        df = pd.read_csv(filepath)
        
        # Handle gene ID column
        if "Unnamed: 0" in df.columns:
            df = df.rename(columns={"Unnamed: 0": "GeneID"})
            df = df.set_index("GeneID")
        elif df.columns[0] in ["GeneID", "gene_id", "Gene"]:
            df = df.set_index(df.columns[0])
        
        # Transpose: genes × samples → samples × genes
        self.raw_counts = df.transpose()
        
        # Store metadata
        self.sample_ids = list(self.raw_counts.index)
        self.gene_ids = list(self.raw_counts.columns)
        self.filepath = filepath
        
        # Performance metrics
        load_time = time.time() - start_time
        memory_mb = self.raw_counts.memory_usage(deep=True).sum() / (1024 * 1024)
        
        self.logger.info(
            f"Loaded in {load_time:.2f}s: {len(self.sample_ids)} samples × "
            f"{len(self.gene_ids)} genes ({memory_mb:.1f} MB)"
        )
        
        return self.raw_counts
    
    def normalize(self, method: str = "cpm") -> pd.DataFrame:
        """
        Normalize expression counts to enable cross-sample comparison.
        
        CPM (Counts Per Million) Formula:
            CPM = (counts / total_counts) × 10^6
        
        This controls for sequencing depth differences between samples.
        
        Args:
            method: Normalization method. Currently supports "cpm".
            
        Returns:
            Normalized expression DataFrame
            
        Raises:
            ValueError: If raw_counts not loaded or invalid method
        """
        if self.raw_counts is None:
            raise ValueError("No data loaded. Call load_data() first.")
        
        self.logger.info(f"Normalizing data using {method.upper()} method")
        
        if method.lower() == "cpm":
            # Calculate library sizes (sum per sample)
            library_sizes = self.raw_counts.sum(axis=1)
            
            # CPM normalization: (count / total) * 1e6
            self.normalized = self.raw_counts.div(library_sizes, axis=0) * 1e6
            
            self.logger.info(
                f"CPM normalization complete. "
                f"Library size range: {library_sizes.min():,.0f} - {library_sizes.max():,.0f}"
            )
        else:
            raise ValueError(f"Unknown normalization method: {method}")
        
        return self.normalized
    
    def filter_genes(
        self,
        min_cpm: float = 1.0,
        min_samples_pct: float = 0.5,
        use_normalized: bool = True
    ) -> List[str]:
        """
        Filter out lowly-expressed genes (noise reduction).
        
        Standard rule: Keep genes with >min_cpm in at least min_samples_pct of samples.
        This removes genes that are essentially not expressed and would add noise
        to downstream analyses.
        
        Args:
            min_cpm: Minimum CPM threshold for a gene to be "expressed"
            min_samples_pct: Minimum fraction of samples where gene must be expressed
            use_normalized: If True, use normalized data; else use raw counts
            
        Returns:
            List of gene IDs that pass the filter
        """
        if use_normalized:
            if self.normalized is None:
                self.logger.warning("Normalized data not found. Running CPM normalization.")
                self.normalize()
            data = self.normalized
        else:
            data = self.raw_counts
        
        n_samples = len(self.sample_ids)
        min_samples = int(n_samples * min_samples_pct)
        
        # Count samples where each gene exceeds threshold
        genes_expressed = (data > min_cpm).sum(axis=0)
        
        # Filter genes
        passing_genes = genes_expressed[genes_expressed >= min_samples].index.tolist()
        
        self.filtered_genes = passing_genes
        n_removed = len(self.gene_ids) - len(passing_genes)
        
        self.logger.info(
            f"Gene filtering: Kept {len(passing_genes):,} / {len(self.gene_ids):,} genes "
            f"(removed {n_removed:,} low-expression genes)"
        )
        self.logger.info(
            f"Filter criteria: >{min_cpm} CPM in ≥{min_samples_pct*100:.0f}% of samples "
            f"(≥{min_samples} samples)"
        )
        
        return self.filtered_genes
    
    def get_filtered_data(self, normalized: bool = True) -> pd.DataFrame:
        """
        Get expression data for filtered genes only.
        
        Args:
            normalized: If True, return normalized data; else raw counts
            
        Returns:
            Filtered expression DataFrame
        """
        if not self.filtered_genes:
            self.logger.warning("No genes filtered yet. Returning all genes.")
            return self.normalized if normalized else self.raw_counts
        
        data = self.normalized if normalized else self.raw_counts
        return data[self.filtered_genes]
    
    def calculate_qc_metrics(self) -> pd.DataFrame:
        """
        Calculate comprehensive quality control metrics for each sample.
        
        Metrics calculated:
        - Library_Size: Total sum of reads per sample
        - Detected_Genes: Count of genes with non-zero expression
        - MT_Percentage: % of reads from mitochondrial genes (MT-*)
        - Log10_Library_Size: Log-transformed library size for visualization
        
        Returns:
            DataFrame with QC metrics for each sample
        """
        if self.raw_counts is None:
            raise ValueError("No data loaded. Call load_data() first.")
        
        self.logger.info("Calculating QC metrics...")
        
        # Library Size: Total reads per sample
        library_sizes = self.raw_counts.sum(axis=1)
        
        # Detected Genes: Non-zero gene count
        detected_genes = (self.raw_counts > 0).sum(axis=1)
        
        # Mitochondrial content
        mt_pct = self.calculate_mt_percentage()
        
        # Create QC DataFrame
        self.qc_metrics = pd.DataFrame({
            "Library_Size": library_sizes,
            "Log10_Library_Size": np.log10(library_sizes + 1),
            "Detected_Genes": detected_genes,
            "MT_Percentage": mt_pct
        })
        
        # Summary statistics
        self.logger.info(f"Average Library Size: {library_sizes.mean():,.0f} reads")
        self.logger.info(f"Average Detected Genes: {detected_genes.mean():,.0f}")
        self.logger.info(f"Average MT%: {mt_pct.mean():.2f}%")
        
        return self.qc_metrics
    
    def calculate_mt_percentage(self) -> pd.Series:
        """
        Calculate the percentage of reads mapping to mitochondrial genes.
        
        High MT content (>20%) often indicates cell membrane rupture during
        sample preparation, leading to cytoplasmic RNA loss while retaining
        mitochondrial RNA. This is a key QC metric for identifying degraded samples.
        
        Returns:
            Series with MT percentage for each sample
        """
        # Find mitochondrial genes (MT- prefix is standard)
        mt_genes = [g for g in self.gene_ids if g.upper().startswith("MT-")]
        
        if not mt_genes:
            self.logger.warning(
                "No mitochondrial genes (MT-*) found. MT% will be 0."
            )
            return pd.Series(0.0, index=self.sample_ids)
        
        self.logger.info(f"Found {len(mt_genes)} mitochondrial genes")
        
        # Calculate MT counts and percentage
        mt_counts = self.raw_counts[mt_genes].sum(axis=1)
        total_counts = self.raw_counts.sum(axis=1)
        mt_percentage = (mt_counts / total_counts) * 100
        
        return mt_percentage
    
    def run_pca(
        self,
        n_components: int = 2,
        use_filtered: bool = True,
        scale: bool = True
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Perform Principal Component Analysis on expression data.
        
        PCA reduces the high-dimensional gene expression data to a small number
        of principal components, allowing visualization of sample relationships.
        It reveals if samples cluster by biological condition or show batch effects.
        
        Args:
            n_components: Number of PCA components to compute
            use_filtered: If True, use filtered genes only
            scale: If True, standardize features before PCA
            
        Returns:
            Tuple of (coordinates array, explained variance ratios)
        """
        # Select data
        if use_filtered and self.filtered_genes:
            data = self.get_filtered_data(normalized=True)
            self.logger.info(f"Running PCA on {len(self.filtered_genes)} filtered genes")
        else:
            if self.normalized is None:
                self.normalize()
            data = self.normalized
            self.logger.info(f"Running PCA on all {len(self.gene_ids)} genes")
        
        # Handle missing values
        data_clean = data.fillna(0)
        
        # Optional scaling
        if scale:
            scaler = StandardScaler()
            data_scaled = scaler.fit_transform(data_clean)
        else:
            data_scaled = data_clean.values
        
        # Run PCA
        pca = PCA(n_components=n_components)
        coordinates = pca.fit_transform(data_scaled)
        
        # Store results
        self.pca_results = {
            "coordinates": coordinates,
            "explained_variance": pca.explained_variance_ratio_,
            "components": pca.components_,
            "n_components": n_components
        }
        
        total_var = sum(pca.explained_variance_ratio_) * 100
        self.logger.info(
            f"PCA complete: {n_components} components explain {total_var:.1f}% of variance"
        )
        for i, var in enumerate(pca.explained_variance_ratio_):
            self.logger.info(f"  PC{i+1}: {var*100:.1f}%")
        
        return coordinates, pca.explained_variance_ratio_
    
    def flag_outliers(
        self,
        lib_size_min: int = 20_000_000,
        lib_size_max: Optional[int] = None,
        mt_threshold: float = 20.0,
        detected_genes_min: Optional[int] = None
    ) -> pd.DataFrame:
        """
        Identify samples that fail QC thresholds.
        
        Args:
            lib_size_min: Minimum acceptable library size
            lib_size_max: Maximum acceptable library size (optional)
            mt_threshold: Maximum acceptable MT percentage
            detected_genes_min: Minimum acceptable detected genes (optional)
            
        Returns:
            DataFrame of failed samples with failure reasons
        """
        if self.qc_metrics is None:
            self.calculate_qc_metrics()
        
        qc = self.qc_metrics.copy()
        qc["Failed"] = False
        qc["Failure_Reasons"] = ""
        
        # Check library size minimum
        lib_low = qc["Library_Size"] < lib_size_min
        qc.loc[lib_low, "Failed"] = True
        qc.loc[lib_low, "Failure_Reasons"] += f"Low library size (<{lib_size_min:,}); "
        
        # Check library size maximum (if specified)
        if lib_size_max:
            lib_high = qc["Library_Size"] > lib_size_max
            qc.loc[lib_high, "Failed"] = True
            qc.loc[lib_high, "Failure_Reasons"] += f"High library size (>{lib_size_max:,}); "
        
        # Check MT percentage
        mt_high = qc["MT_Percentage"] > mt_threshold
        qc.loc[mt_high, "Failed"] = True
        qc.loc[mt_high, "Failure_Reasons"] += f"High MT% (>{mt_threshold}%); "
        
        # Check detected genes (if specified)
        if detected_genes_min:
            genes_low = qc["Detected_Genes"] < detected_genes_min
            qc.loc[genes_low, "Failed"] = True
            qc.loc[genes_low, "Failure_Reasons"] += f"Low gene detection (<{detected_genes_min:,}); "
        
        # Get failed samples
        failed = qc[qc["Failed"]][["Library_Size", "MT_Percentage", "Detected_Genes", "Failure_Reasons"]]
        
        n_failed = len(failed)
        n_total = len(qc)
        
        self.logger.info(
            f"Outlier detection: {n_failed}/{n_total} samples failed QC "
            f"({n_failed/n_total*100:.1f}%)"
        )
        
        if n_failed > 0:
            self.logger.warning(f"Failed samples: {list(failed.index)}")
        
        return failed
    
    def save_failed_samples(
        self,
        output_path: str,
        lib_size_min: int = 20_000_000,
        mt_threshold: float = 20.0
    ) -> str:
        """
        Export failed samples to CSV file.
        
        Args:
            output_path: Path to save the CSV file
            lib_size_min: Minimum library size threshold
            mt_threshold: Maximum MT percentage threshold
            
        Returns:
            Path to saved file
        """
        failed = self.flag_outliers(
            lib_size_min=lib_size_min,
            mt_threshold=mt_threshold
        )
        
        if len(failed) > 0:
            failed.to_csv(output_path)
            self.logger.info(f"Saved {len(failed)} failed samples to: {output_path}")
        else:
            # Create empty file with header
            pd.DataFrame(
                columns=["Library_Size", "MT_Percentage", "Detected_Genes", "Failure_Reasons"]
            ).to_csv(output_path)
            self.logger.info(f"No failed samples. Empty file saved to: {output_path}")
        
        return output_path
    
    def get_summary(self) -> Dict[str, Any]:
        """
        Get a summary of the dataset and analysis state.
        
        Returns:
            Dictionary with dataset summary information
        """
        summary = {
            "filepath": self.filepath,
            "n_samples": len(self.sample_ids),
            "n_genes_total": len(self.gene_ids),
            "n_genes_filtered": len(self.filtered_genes) if self.filtered_genes else None,
            "normalized": self.normalized is not None,
            "qc_calculated": self.qc_metrics is not None,
            "pca_run": self.pca_results is not None
        }
        
        if self.qc_metrics is not None:
            summary["avg_library_size"] = self.qc_metrics["Library_Size"].mean()
            summary["avg_detected_genes"] = self.qc_metrics["Detected_Genes"].mean()
            summary["avg_mt_pct"] = self.qc_metrics["MT_Percentage"].mean()
        
        return summary
    
    def __repr__(self) -> str:
        """String representation of the dataset."""
        status = []
        if self.normalized is not None:
            status.append("normalized")
        if self.filtered_genes:
            status.append(f"filtered({len(self.filtered_genes)} genes)")
        if self.qc_metrics is not None:
            status.append("qc_computed")
        
        status_str = ", ".join(status) if status else "raw"
        
        return (
            f"RNASeqDataset({len(self.sample_ids)} samples × {len(self.gene_ids)} genes, "
            f"status=[{status_str}])"
        )
