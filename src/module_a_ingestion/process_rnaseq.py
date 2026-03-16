import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def process_differential_expression(filepath, pval_cutoff=0.05):
    """Ingests RNA-seq data and maps it to the 0.01 - 100 solver bounds."""
    df = pd.read_csv(filepath)
    
    # 1. Biological Filtering
    df['is_significant'] = df['padj'] <= pval_cutoff
    
    # 2. Mathematical Mapping: log2(FC) -> Absolute FC
    df['mp_biopath_input'] = (2 ** df['log2FoldChange']).clip(0.01, 100.0)
    
    return df

def plot_mapping_results(df, pval_cutoff=0.05):
    """Generates a volcano plot highlighting the MP-BioPath mapping."""
    plt.figure(figsize=(10, 6))
    
    # Plot non-significant noise in gray
    noise = df[~df['is_significant']]
    plt.scatter(noise['log2FoldChange'], -np.log10(noise['padj']), c='gray', alpha=0.5, label='Insignificant Noise')
    
    # Plot significant hits in blue
    hits = df[df['is_significant']]
    plt.scatter(hits['log2FoldChange'], -np.log10(hits['padj']), c='blue', alpha=0.7, label='MP-BioPath Inputs')
    
    plt.axhline(-np.log10(pval_cutoff), color='red', linestyle='--', label='p=0.05 Threshold')
    plt.title("Genomic Data Mapping: log2(FC) to MP-BioPath Boundaries")
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10(adj p-value)")
    plt.legend()
    plt.savefig("../../docs/mapping_volcano.png")
    print("Visualization saved to docs/mapping_volcano.png")

if __name__ == "__main__":
    data = process_differential_expression("../../data/mock_rnaseq.csv")
    plot_mapping_results(data)
    print(data[['Gene', 'mp_biopath_input']])