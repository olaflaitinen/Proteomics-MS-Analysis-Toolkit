# python/src/enrichment.py

import gseapy as gp
import pandas as pd

def run_gsea(de_results, gene_sets='GO_Biological_Process_2021', organism='Human'):
    """
    Performs Gene Set Enrichment Analysis (GSEA) using gseapy.

    Args:
        de_results (pd.DataFrame): DataFrame from differential analysis,
                                   must have log2FC as a column and gene names as index.
        gene_sets (str or list): Name of the gene set library to use (e.g., from Enrichr).
        organism (str): Organism for the gene sets.

    Returns:
        pd.DataFrame: The GSEA results.
    """
    # Create a ranked list of genes based on log2FC
    ranked_list = de_results.sort_values('log2FC', ascending=False)[['log2FC']]
    
    print(f"Running GSEA with {len(ranked_list)} ranked genes...")
    
    # Run pre-ranked GSEA
    pre_res = gp.prerank(
        rnk=ranked_list,
        gene_sets=gene_sets,
        organism=organism,
        seed=42
    )
    
    return pre_res.res2d
