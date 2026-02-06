import pandas as pd
import numpy as np

def process_expression_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Cleans the raw expression matrix:
    1. Sets Gene IDs as the index.
    2. Transposes (Rows=Patients, Cols=Genes).
    """
    # The first column is "Unnamed: 0" (Gene IDs). Let's fix that.
    if "Unnamed: 0" in df.columns:
        df = df.rename(columns={"Unnamed: 0": "GeneID"})
        df = df.set_index("GeneID")
    
    # Transpose: Flip the matrix 90 degrees
    # Now: Rows = Patients, Cols = Genes
    df_transposed = df.transpose()
    
    print(f"--- Data Transposed ---")
    print(f"New Dimensions: {df_transposed.shape[0]} Patients x {df_transposed.shape[1]} Genes")
    
    return df_transposed

def calculate_qc_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates biological quality metrics for each patient.
    """
    print("\n--- Calculating QC Metrics ---")
    
    # 1. Library Size: Total sum of reads per patient
    # Biologically: Did we sequence this patient deeply enough?
    library_sizes = df.sum(axis=1)
    
    # 2. Detected Genes: Count of genes with non-zero counts
    # Biologically: Is the library complex or mostly empty noise?
    detected_genes = (df > 0).sum(axis=1)
    
    # Create a clean QC table
    qc_df = pd.DataFrame({
        "Library_Size": library_sizes,
        "Detected_Genes": detected_genes
    })
    
    # Basic Stats
    mean_lib_size = qc_df['Library_Size'].mean()
    print(f"Average Library Size: {mean_lib_size:,.0f} reads")
    
    return qc_df

if __name__ == "__main__":
    from src.loader import load_real_data
    
    # 1. Load (Using the path to the big file)
    FILE_PATH = "data/TCGA_LIHC_Gene_Expression.csv"
    df = load_real_data(FILE_PATH)
    
    # 2. Process (Transpose)
    df_clean = process_expression_data(df)
    
    # 3. Analyze (QC)
    qc_stats = calculate_qc_metrics(df_clean)
    
    # 4. Preview Bad Samples (Optional: Show samples with low counts)
    print("\nTop 5 Samples with Lowest Library Size:")
    print(qc_stats.sort_values("Library_Size").head(5))