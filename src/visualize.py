import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

def plot_library_sizes(qc_df: pd.DataFrame, output_dir: str):
    """
    Generates a histogram of library sizes.
    Answers: "Do all patients have enough data?"
    """
    plt.figure(figsize=(10, 6))
    
    # Create the Histogram
    sns.histplot(qc_df['Library_Size'], bins=30, kde=True, color="teal")
    
    # Add Threshold Line (e.g., 20 Million Reads)
    plt.axvline(20_000_000, color='red', linestyle='--', label='Min Recommended (20M)')
    
    plt.title("Distribution of Library Sizes (Total Reads per Patient)")
    plt.xlabel("Library Size (Reads)")
    plt.ylabel("Count of Patients")
    plt.legend()
    
    # Save
    save_path = os.path.join(output_dir, "qc_library_sizes.png")
    plt.savefig(save_path)
    print(f"[Saved] Library Size Plot: {save_path}")
    plt.close()

def plot_detected_genes(qc_df: pd.DataFrame, output_dir: str):
    """
    Generates a boxplot of detected genes.
    Answers: "Are any samples failing to detect genes?"
    """
    plt.figure(figsize=(10, 6))
    
    # Create Boxplot
    sns.boxplot(y=qc_df['Detected_Genes'], color="cornflowerblue")
    
    plt.title("Distribution of Detected Genes per Patient")
    plt.ylabel("Number of Genes Detected (>0 counts)")
    
    # Save
    save_path = os.path.join(output_dir, "qc_detected_genes.png")
    plt.savefig(save_path)
    print(f"[Saved] Detected Genes Plot: {save_path}")
    plt.close()

if __name__ == "__main__":
    from src.loader import load_real_data
    from src.qc import process_expression_data, calculate_qc_metrics
    
    # 1. Pipeline: Load -> Process -> QC
    FILE_PATH = "data/TCGA_LIHC_Gene_Expression.csv"
    OUTPUT_DIR = "output"
    
    # Ensure output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    df = load_real_data(FILE_PATH)
    df_clean = process_expression_data(df)
    qc_stats = calculate_qc_metrics(df_clean)
    
    # 2. Visualize
    print("\n--- Generating Dashboard Plots ---")
    plot_library_sizes(qc_stats, OUTPUT_DIR)
    plot_detected_genes(qc_stats, OUTPUT_DIR)