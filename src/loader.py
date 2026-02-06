import pandas as pd
import time
import os

def load_real_data(filepath: str) -> pd.DataFrame:
    """
    Loads a large dataset and profiles memory usage.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Error: The file '{filepath}' was not found.")

    print(f"--- Loading Real Data: {filepath} ---")
    
    # Start timer for performance measurement
    start_time = time.time()
    
    # Load the CSV
    df = pd.read_csv(filepath)
    
    end_time = time.time()
    
    # 1. Measure Scale
    print(f"Loaded in {end_time - start_time:.4f} seconds")
    print(f"Total Samples (Rows): {df.shape[0]}")
    print(f"Total Features (Columns): {df.shape[1]}")
    
    # 2. Memory Usage Profiling
    memory_bytes = df.memory_usage(deep=True).sum()
    memory_mb = memory_bytes / (1024 * 1024)
    print(f"Memory Usage: {memory_mb:.2f} MB")
    
    return df

if __name__ == "__main__":
    # Point to data
    FILE_PATH = "data/TCGA_LIHC_Gene_Expression.csv"
    try:
        df = load_real_data(FILE_PATH)
        print("\nFirst 5 rows (Preview only):")
        print(df.head())
        
    except Exception as e:
        print(f"An error occurred: {e}")