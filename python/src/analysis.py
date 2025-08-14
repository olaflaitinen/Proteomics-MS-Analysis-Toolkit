# python/src/analysis.py

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

def run_differential_analysis(processed_data, metadata_path, control_group='control', treat_group='treatment'):
    """
    Performs differential abundance analysis using Welch's t-test.

    Args:
        processed_data (pd.DataFrame): Normalized and imputed data matrix.
        metadata_path (str): Path to the metadata CSV file.
        control_group (str): Name of the control condition.
        treat_group (str): Name of the treatment condition.

    Returns:
        pd.DataFrame: A dataframe with log2FC, p-value, and adjusted p-value.
    """
    meta = pd.read_csv(metadata_path)
    
    # Get sample names for each group
    control_samples = meta[meta['condition'] == control_group]['sample'].tolist()
    treat_samples = meta[meta['condition'] == treat_group]['sample'].tolist()
    
    # Format column names to match data
    control_cols = [f'Intensity {s}' for s in control_samples]
    treat_cols = [f'Intensity {s}' for s in treat_samples]
    
    # Calculate log2 Fold Change
    control_mean = processed_data[control_cols].mean(axis=1)
    treat_mean = processed_data[treat_cols].mean(axis=1)
    log2fc = treat_mean - control_mean
    
    # Perform Welch's t-test for each protein
    p_values = ttest_ind(
        processed_data[treat_cols], 
        processed_data[control_cols], 
        axis=1, 
        equal_var=False, # Welch's t-test
        nan_policy='omit'
    ).pvalue
    
    # Create results dataframe
    results_df = pd.DataFrame({
        'log2FC': log2fc,
        'p_value': p_values
    }, index=processed_data.index)
    
    # Simple Benjamini-Hochberg adjustment for multiple testing
    results_df = results_df.sort_values('p_value')
    results_df['adj_p_value'] = results_df['p_value'] * len(results_df) / np.arange(1, len(results_df) + 1)
    results_df['adj_p_value'] = np.minimum(results_df['adj_p_value'], 1.0) # Cap at 1.0
    
    return results_df.sort_index()
