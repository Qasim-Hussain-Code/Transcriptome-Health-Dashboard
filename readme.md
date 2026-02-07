# Transcriptome Health Dashboard

A Bioinformatics QC pipeline for assessing RNA-Seq data quality. This tool profiles library sizes, detects gene counts, and generates visual dashboards to identify failed samples.

## Features
* **Big Data Loader:** Memory-efficient loading of large expression matrices (60k+ genes).
* **Biological QC:** Calculates Library Size and Gene Detection Rates.
* **Visualization:** Generates publication-ready distributions.

## Usage

```bash
# 1. Load Data & Run QC
python -m src.qc

# 2. Generate Dashboard Plots
python -m src.visualize