# Transcriptome Health Dashboard v2.0

A **professional-grade** RNA-Seq quality control and analysis pipeline with biologically-rigorous normalization, interactive visualizations, and comprehensive sample QC metrics.

![Python](https://img.shields.io/badge/Python-3.9+-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![Tests](https://img.shields.io/badge/Tests-Passing-brightgreen.svg)
=======
## Project Overview
* **Dataset:** TCGA-LIHC (Liver Hepatocellular Carcinoma)
* **Scale:** 60,660 Genes x 425 Patients
* **Goal:** Diagnose sequencing depth and library complexity before downstream AI analysis.

## Key Results

### 1. Library Size Distribution
*This histogram proves that all patients exceed the minimum threshold of 20 Million reads (Red Line). The "Bell Curve" shape indicates consistent sequencing depth across the cohort.*

![Library Size Distribution](output/qc_library_sizes.png)


---

## Overview

This pipeline performs comprehensive quality control analysis for bulk RNA-Seq datasets, implementing industry-standard bioinformatics practices:

| Feature | Description |
|---------|-------------|
| **CPM Normalization** | Counts Per Million normalization for cross-sample comparability |
| **Gene Filtering** | Remove lowly-expressed genes (noise reduction) |
| **PCA Analysis** | Visualize sample structure and detect batch effects |
| **Mitochondrial QC** | Flag degraded samples by MT-gene content |
| **Outlier Detection** | Automatic flagging of failed samples |
| **Interactive Dashboards** | Hover-enabled Plotly HTML reports |

### Dataset

* **Source:** TCGA-LIHC (Liver Hepatocellular Carcinoma)
* **Scale:** 60,660 Genes × 425 Patients
* **Goal:** Diagnose sequencing depth and library complexity before downstream analysis

---

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/Qasim-Hussain/Transcriptome-Health-Dashboard.git
cd Transcriptome-Health-Dashboard

<<<<<<< HEAD
# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Usage

```bash
# Run the complete QC pipeline
python -m src.main --input data/TCGA_LIHC_Gene_Expression.csv --output results/

# With custom thresholds
python -m src.main \
  --input data/TCGA_LIHC_Gene_Expression.csv \
  --output results/ \
  --lib-threshold 25000000 \
  --mt-threshold 15
```

### Python API

```python
from src.dataset import RNASeqDataset

# Load and analyze
dataset = RNASeqDataset("data/TCGA_LIHC_Gene_Expression.csv")
dataset.normalize(method="cpm")
dataset.filter_genes(min_cpm=1.0, min_samples_pct=0.5)
dataset.calculate_qc_metrics()
dataset.run_pca(n_components=2)

# Export failed samples
dataset.save_failed_samples("output/failed_samples.csv")

# Generate interactive dashboard
from src.visualize import create_interactive_dashboard
create_interactive_dashboard(dataset, output_dir="output/")
```

---

## QC Metrics Explained

### 1. Library Size Distribution
**Question:** Do all samples have sufficient sequencing depth?

> The library size is the total count of reads per sample. Samples below 20M reads may have insufficient coverage for reliable gene expression quantification.

### 2. Gene Detection Complexity
**Question:** Are libraries diverse or dominated by a few genes?

> Samples should detect >15,000 genes. Low detection may indicate RNA degradation or library preparation issues.

### 3. Mitochondrial Content (MT%)
**Question:** Are samples degraded?

> High mitochondrial content (>20%) often indicates cell membrane rupture during sample preparation, leading to cytoplasmic RNA loss while retaining mitochondrial RNA. This is a key QC metric.

**Formula:** `MT% = (∑ MT-gene counts / Total counts) × 100`

### 4. PCA - Sample Structure
**Question:** Do samples cluster by biological condition?

> PCA reduces 60,000+ genes to 2 dimensions. Samples should cluster by disease status, not by technical factors (batch effect).

---

## 📐 Scientific Methods

### CPM Normalization

Counts Per Million (CPM) normalization controls for sequencing depth:

$$\text{CPM} = \frac{\text{counts}}{\text{total counts}} \times 10^6$$

**Why?** You cannot compare Gene X in Patient A vs. Patient B if Patient A has 50M reads and Patient B has 20M. CPM makes samples comparable.

### Gene Filtering

Low-expression genes add noise without biological signal:

**Rule:** Keep genes with `>1 CPM` in `≥50%` of samples

This typically reduces 60,000 genes to ~15,000-20,000 informative genes.

---

## Project Structure

```
Transcriptome-Health-Dashboard/
├── src/
│   ├── __init__.py
│   ├── dataset.py      # Core RNASeqDataset class
│   ├── main.py         # CLI entry point
│   ├── visualize.py    # Interactive Plotly visualizations
│   ├── qc.py           # Legacy QC functions
│   └── loader.py       # Legacy data loading
├── tests/
│   ├── conftest.py     # Pytest fixtures
│   ├── test_normalization.py
│   ├── test_filtering.py
│   ├── test_qc_metrics.py
│   └── test_pca.py
├── data/
│   └── TCGA_LIHC_Gene_Expression.csv
├── output/
│   ├── dashboard.html           # Interactive QC dashboard
│   ├── failed_samples.csv       # Samples failing QC
│   ├── qc_metrics.csv           # Per-sample metrics
│   └── pca_coordinates.csv      # PCA results
├── requirements.txt
├── pyproject.toml
└── README.md
```

---

## CLI Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--input`, `-i` | *Required* | Path to expression CSV |
| `--output`, `-o` | *Required* | Output directory |
| `--lib-threshold` | 20,000,000 | Minimum library size |
| `--mt-threshold` | 20.0 | Maximum MT percentage |
| `--min-cpm` | 1.0 | CPM threshold for gene filtering |
| `--min-samples-pct` | 0.5 | Fraction of samples for filtering |
| `--skip-pca` | False | Skip PCA analysis |
| `--log-level` | INFO | Logging verbosity |

---

## Testing

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=src --cov-report=html

# Run specific test file
pytest tests/test_normalization.py -v
```

---

## Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| pandas | 2.1.4 | Data manipulation |
| numpy | 1.26.2 | Numerical operations |
| scikit-learn | 1.3.2 | PCA, StandardScaler |
| plotly | 5.18.0 | Interactive visualizations |
| matplotlib | 3.8.2 | Static plots |
| seaborn | 0.13.0 | Statistical visualizations |
| pytest | 7.4.3 | Unit testing |

---

## Output Files

| File | Description |
|------|-------------|
| `dashboard.html` | Interactive QC dashboard with all plots |
| `qc_metrics.csv` | Library size, detected genes, MT% per sample |
| `failed_samples.csv` | Samples failing QC thresholds with reasons |
| `pca_coordinates.csv` | PC1, PC2 coordinates for each sample |
| `pipeline.log` | Detailed execution log |
| `*.html` | Individual interactive plots |
| `*.png` | Static plots (backward compatibility) |

---

## Key Results

### Library Size Distribution
*All patients exceed the minimum threshold of 20M reads. The bell-curve shape indicates consistent sequencing depth.*

![Library Size](output/qc_library_sizes.png)

### Gene Detection Complexity
*No samples show extreme dropout, confirming high-quality libraries.*

![Detected Genes](output/qc_detected_genes.png)

---

## References

1. **CPM Normalization**: Robinson MD, Oshlack A. A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biology* (2010).

2. **Gene Filtering**: Chen Y, Lun AT, Smyth GK. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. *F1000Research* (2016).

3. **MT Content as QC**: Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems Biology* (2019).

---

## License

MIT License - See [LICENSE](LICENSE) for details.

---
