import pandas as pd
import numpy as np

def process_differential_expression(filepath, pval_cutoff=0.05):
    """
    Ingests RNA-seq data and maps it to the 0.01 - 100 bounds required 
    by the MP-BioPath optimization solver.
    """
    print(f"Loading genomic data from {filepath}...")
    df = pd.read_csv(filepath)
    
    # 1. Biological Filtering: Remove statistically insignificant noise
    # We only want to feed actual perturbations to the solver.
    initial_count = len(df)
    df = df[df['padj'] <= pval_cutoff].copy()
    filtered_count = len(df)
    print(f"Filtered out {initial_count - filtered_count} insignificant genes (padj > {pval_cutoff}).")

    # 2. Mathematical Mapping: log2(FC) to Absolute FC
    # MP-BioPath baseline is 1.0. A log2FC of 0 becomes 1.0.
    df['absolute_fc'] = 2 ** df['log2FoldChange']

    # 3. Solver Constraints: Clip values to [0.01, 100]
    # MP-BioPath equations fail if bounds are exceeded.
    df['mp_biopath_input'] = df['absolute_fc'].clip(lower=0.01, upper=100.0)

    return df[['Gene', 'log2FoldChange', 'mp_biopath_input']]

if __name__ == "__main__":
    # Execute the pipeline on our test data
    input_file = "../../data/mock_rnaseq.csv"
    output_data = process_differential_expression(input_file)
    
    print("\n--- Final MP-BioPath Solver Inputs ---")
    print(output_data.to_string(index=False))