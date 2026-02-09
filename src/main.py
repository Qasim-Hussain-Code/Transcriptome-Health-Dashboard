"""
Command Line Interface for RNA-Seq Quality Control Pipeline

This module provides a professional CLI entry point for the Transcriptome
Health Dashboard, allowing users to run the complete QC pipeline with
configurable parameters.

Usage:
    python -m src.main --input data/expression.csv --output results/
    python -m src.main --input data/TCGA.csv --output output/ --lib-threshold 25000000 --mt-threshold 15

Arguments:
    --input: Path to expression CSV file (genes × samples, first col = GeneID)
    --output: Directory for output files (created if doesn't exist)
    --lib-threshold: Minimum library size threshold (default: 20,000,000)
    --mt-threshold: Maximum mitochondrial % threshold (default: 20.0)
    --min-cpm: Minimum CPM for gene filtering (default: 1.0)
    --min-samples-pct: Minimum % of samples for gene filtering (default: 0.5)
    --log-level: Logging verbosity (DEBUG, INFO, WARNING, ERROR)
    --skip-pca: Skip PCA analysis
    --skip-interactive: Generate only static PNG plots (no HTML dashboard)
"""

import argparse
import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

# Add parent to path for module imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.dataset import RNASeqDataset
from src.visualize import (
    create_interactive_dashboard,
    save_individual_plots,
    plot_library_sizes_interactive,
    plot_detected_genes_interactive,
    plot_mt_content_interactive,
    plot_pca_interactive
)


def setup_logging(
    log_level: str = "INFO",
    log_file: Optional[str] = None
) -> logging.Logger:
    """
    Configure professional logging with both console and file output.
    
    Args:
        log_level: Logging verbosity level
        log_file: Optional path to log file
        
    Returns:
        Configured logger instance
    """
    # Create logger
    logger = logging.getLogger("RNASeqPipeline")
    logger.setLevel(getattr(logging, log_level.upper()))
    
    # Clear existing handlers
    logger.handlers = []
    
    # Console handler with colored output
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG)
    console_format = logging.Formatter(
        '%(asctime)s │ %(levelname)-8s │ %(message)s',
        datefmt='%H:%M:%S'
    )
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)
    
    # File handler (if specified)
    if log_file:
        file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_format = logging.Formatter(
            '%(asctime)s | %(levelname)-8s | %(name)s | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(file_format)
        logger.addHandler(file_handler)
    
    return logger


def create_argument_parser() -> argparse.ArgumentParser:
    """Create the CLI argument parser with all options."""
    parser = argparse.ArgumentParser(
        prog="transcriptome-qc",
        description="""
╔══════════════════════════════════════════════════════════════════╗
║           TRANSCRIPTOME HEALTH DASHBOARD v2.0                    ║
║    Professional RNA-Seq Quality Control & Analysis Pipeline      ║
╚══════════════════════════════════════════════════════════════════╝

This pipeline performs comprehensive QC analysis including:
  • CPM normalization for cross-sample comparability
  • Gene filtering to remove low-expression noise
  • Library size and complexity assessment
  • Mitochondrial content detection (cell degradation marker)
  • PCA for sample structure visualization
  • Automated outlier detection and flagging
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Basic usage:
    python -m src.main --input data/expression.csv --output results/
    
  With custom thresholds:
    python -m src.main -i data/TCGA.csv -o output/ --lib-threshold 25000000 --mt-threshold 15
    
  Debug mode with all logging:
    python -m src.main -i data/data.csv -o out/ --log-level DEBUG
        """
    )
    
    # Required arguments
    required = parser.add_argument_group("Required Arguments")
    required.add_argument(
        "-i", "--input",
        type=str,
        required=True,
        metavar="FILE",
        help="Path to expression CSV file (genes × samples format)"
    )
    required.add_argument(
        "-o", "--output",
        type=str,
        required=True,
        metavar="DIR",
        help="Output directory for results (created if doesn't exist)"
    )
    
    # QC threshold arguments
    thresholds = parser.add_argument_group("QC Thresholds")
    thresholds.add_argument(
        "--lib-threshold",
        type=int,
        default=20_000_000,
        metavar="N",
        help="Minimum library size threshold (default: 20,000,000)"
    )
    thresholds.add_argument(
        "--mt-threshold",
        type=float,
        default=20.0,
        metavar="PCT",
        help="Maximum mitochondrial %% threshold (default: 20.0)"
    )
    
    # Gene filtering arguments
    filtering = parser.add_argument_group("Gene Filtering")
    filtering.add_argument(
        "--min-cpm",
        type=float,
        default=1.0,
        metavar="CPM",
        help="Minimum CPM for gene to be considered expressed (default: 1.0)"
    )
    filtering.add_argument(
        "--min-samples-pct",
        type=float,
        default=0.5,
        metavar="PCT",
        help="Minimum fraction of samples where gene must be expressed (default: 0.5)"
    )
    
    # Analysis options
    options = parser.add_argument_group("Analysis Options")
    options.add_argument(
        "--skip-pca",
        action="store_true",
        help="Skip PCA analysis"
    )
    options.add_argument(
        "--skip-interactive",
        action="store_true",
        help="Generate only static PNG plots (no HTML dashboard)"
    )
    options.add_argument(
        "--pca-components",
        type=int,
        default=2,
        metavar="N",
        help="Number of PCA components to compute (default: 2)"
    )
    
    # Logging
    logging_group = parser.add_argument_group("Logging")
    logging_group.add_argument(
        "--log-level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging verbosity (default: INFO)"
    )
    
    return parser


def print_banner(logger: logging.Logger):
    """Print the pipeline banner."""
    banner = """
    ╔══════════════════════════════════════════════════════════════════╗
    ║           TRANSCRIPTOME HEALTH DASHBOARD v2.0                    ║
    ║    Professional RNA-Seq Quality Control & Analysis Pipeline      ║
    ╚══════════════════════════════════════════════════════════════════╝
    """
    logger.info(banner)


def run_pipeline(args: argparse.Namespace, logger: logging.Logger) -> int:
    """
    Execute the complete QC pipeline.
    
    Args:
        args: Parsed command line arguments
        logger: Logger instance
        
    Returns:
        Exit code (0 = success, 1 = error)
    """
    start_time = datetime.now()
    
    try:
        # Create output directory
        os.makedirs(args.output, exist_ok=True)
        logger.info(f"Output directory: {args.output}")
        
        # ═══════════════════════════════════════════════════════════════
        # STEP 1: Load Data
        # ═══════════════════════════════════════════════════════════════
        logger.info("═" * 60)
        logger.info("STEP 1: Loading expression data")
        logger.info("═" * 60)
        
        dataset = RNASeqDataset(filepath=args.input, logger=logger)
        
        # ═══════════════════════════════════════════════════════════════
        # STEP 2: Normalize (CPM)
        # ═══════════════════════════════════════════════════════════════
        logger.info("═" * 60)
        logger.info("STEP 2: CPM Normalization")
        logger.info("═" * 60)
        
        dataset.normalize(method="cpm")
        
        # ═══════════════════════════════════════════════════════════════
        # STEP 3: Gene Filtering
        # ═══════════════════════════════════════════════════════════════
        logger.info("═" * 60)
        logger.info("STEP 3: Gene Filtering (Noise Reduction)")
        logger.info("═" * 60)
        
        dataset.filter_genes(
            min_cpm=args.min_cpm,
            min_samples_pct=args.min_samples_pct
        )
        
        # ═══════════════════════════════════════════════════════════════
        # STEP 4: Calculate QC Metrics
        # ═══════════════════════════════════════════════════════════════
        logger.info("═" * 60)
        logger.info("STEP 4: Calculating QC Metrics")
        logger.info("═" * 60)
        
        dataset.calculate_qc_metrics()
        
        # Save QC metrics to CSV
        qc_path = os.path.join(args.output, "qc_metrics.csv")
        dataset.qc_metrics.to_csv(qc_path)
        logger.info(f"Saved QC metrics: {qc_path}")
        
        # ═══════════════════════════════════════════════════════════════
        # STEP 5: PCA Analysis (Optional)
        # ═══════════════════════════════════════════════════════════════
        if not args.skip_pca:
            logger.info("═" * 60)
            logger.info("STEP 5: PCA Analysis")
            logger.info("═" * 60)
            
            dataset.run_pca(n_components=args.pca_components)
            
            # Save PCA coordinates
            pca_df = dataset.qc_metrics.copy()
            for i in range(args.pca_components):
                pca_df[f"PC{i+1}"] = dataset.pca_results["coordinates"][:, i]
            pca_path = os.path.join(args.output, "pca_coordinates.csv")
            pca_df.to_csv(pca_path)
            logger.info(f"Saved PCA coordinates: {pca_path}")
        
        # ═══════════════════════════════════════════════════════════════
        # STEP 6: Outlier Detection
        # ═══════════════════════════════════════════════════════════════
        logger.info("═" * 60)
        logger.info("STEP 6: Outlier Detection")
        logger.info("═" * 60)
        
        failed_path = os.path.join(args.output, "failed_samples.csv")
        dataset.save_failed_samples(
            output_path=failed_path,
            lib_size_min=args.lib_threshold,
            mt_threshold=args.mt_threshold
        )
        
        # ═══════════════════════════════════════════════════════════════
        # STEP 7: Generate Visualizations
        # ═══════════════════════════════════════════════════════════════
        logger.info("═" * 60)
        logger.info("STEP 7: Generating Visualizations")
        logger.info("═" * 60)
        
        if not args.skip_interactive:
            # Create interactive HTML dashboard
            dashboard_path = create_interactive_dashboard(
                dataset=dataset,
                output_dir=args.output,
                lib_threshold=args.lib_threshold,
                mt_threshold=args.mt_threshold
            )
            logger.info(f"Interactive dashboard: {dashboard_path}")
            
            # Save individual plots as HTML
            plot_paths = save_individual_plots(
                dataset=dataset,
                output_dir=args.output,
                lib_threshold=args.lib_threshold,
                mt_threshold=args.mt_threshold
            )
            for name, path in plot_paths.items():
                logger.info(f"  - {name}: {path}")
        
        # Also save static PNGs for compatibility
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            
            # Static library size histogram
            plt.figure(figsize=(10, 6))
            sns.histplot(dataset.qc_metrics['Library_Size'], bins=30, kde=True, color="teal")
            plt.axvline(args.lib_threshold, color='red', linestyle='--', 
                       label=f'Min Threshold ({args.lib_threshold/1e6:.0f}M)')
            plt.title("Distribution of Library Sizes")
            plt.xlabel("Library Size (Reads)")
            plt.ylabel("Count")
            plt.legend()
            plt.savefig(os.path.join(args.output, "qc_library_sizes.png"), dpi=150, bbox_inches='tight')
            plt.close()
            
            # Static detected genes boxplot
            plt.figure(figsize=(10, 6))
            sns.boxplot(y=dataset.qc_metrics['Detected_Genes'], color="cornflowerblue")
            plt.title("Distribution of Detected Genes per Sample")
            plt.ylabel("Number of Genes Detected")
            plt.savefig(os.path.join(args.output, "qc_detected_genes.png"), dpi=150, bbox_inches='tight')
            plt.close()
            
            logger.info("Static PNG plots saved for backward compatibility")
            
        except ImportError:
            logger.warning("matplotlib/seaborn not available for static plots")
        
        # ═══════════════════════════════════════════════════════════════
        # COMPLETE
        # ═══════════════════════════════════════════════════════════════
        elapsed = datetime.now() - start_time
        
        logger.info("═" * 60)
        logger.info("PIPELINE COMPLETE")
        logger.info("═" * 60)
        logger.info(f"Total time: {elapsed.total_seconds():.2f} seconds")
        logger.info(f"Samples analyzed: {len(dataset.sample_ids)}")
        logger.info(f"Genes (original): {len(dataset.gene_ids)}")
        logger.info(f"Genes (filtered): {len(dataset.filtered_genes)}")
        logger.info(f"Output directory: {args.output}")
        
        # Summary of outputs
        logger.info("\nOutput Files:")
        for f in sorted(os.listdir(args.output)):
            path = os.path.join(args.output, f)
            if os.path.isfile(path):
                size_kb = os.path.getsize(path) / 1024
                logger.info(f"  └── {f} ({size_kb:.1f} KB)")
        
        return 0
        
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        return 1
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        return 1


def main():
    """Main entry point for CLI."""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    # Setup logging
    log_file = os.path.join(args.output, "pipeline.log") if args.output else None
    os.makedirs(args.output, exist_ok=True)
    logger = setup_logging(log_level=args.log_level, log_file=log_file)
    
    # Print banner
    print_banner(logger)
    
    # Log configuration
    logger.info("Configuration:")
    logger.info(f"  Input file: {args.input}")
    logger.info(f"  Output directory: {args.output}")
    logger.info(f"  Library size threshold: {args.lib_threshold:,}")
    logger.info(f"  MT threshold: {args.mt_threshold}%")
    logger.info(f"  Gene filter: >{args.min_cpm} CPM in ≥{args.min_samples_pct*100:.0f}% samples")
    logger.info(f"  Skip PCA: {args.skip_pca}")
    logger.info(f"  Skip interactive: {args.skip_interactive}")
    
    # Run pipeline
    exit_code = run_pipeline(args, logger)
    
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
