# python/src/plotting.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def plot_pca(processed_data, metadata_path):
    """
    Performs PCA and plots the first two principal components.
    """
    meta = pd.read_csv(metadata_path)
    
    # Transpose data so that samples are rows and proteins are columns
    X = processed_data.T
    
    # Scale data
    X_scaled = StandardScaler().fit_transform(X)
    
    # Perform PCA
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(X_scaled)
    pc_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
    
    # Add metadata for plotting
    pc_df['sample'] = X.index
    pc_df = pc_df.merge(meta.rename(columns={'sample': 'mq_sample_name'}), left_on='sample', right_on=pc_df.columns[0], how='left')

    # Plot
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        x='PC1', 
        y='PC2', 
        hue='condition', 
        data=pc_df, 
        s=100,
        alpha=0.8
    )
    plt.title('PCA of Proteomics Data', fontsize=16)
    plt.xlabel(f'Principal Component 1 ({pca.explained_variance_ratio_[0]:.1%})', fontsize=12)
    plt.ylabel(f'Principal Component 2 ({pca.explained_variance_ratio_[1]:.1%})', fontsize=12)
    plt.legend(title='Condition')
    plt.grid(True)
    plt.show()

def plot_volcano(de_results, p_threshold=0.05, fc_threshold=1.0):
    """
    Creates a volcano plot from differential analysis results.
    """
    # Create a column for significance
    de_results['significance'] = 'Not Significant'
    de_results.loc[(de_results['adj_p_value'] < p_threshold) & (de_results['log2FC'] > fc_threshold), 'significance'] = 'Upregulated'
    de_results.loc[(de_results['adj_p_value'] < p_threshold) & (de_results['log2FC'] < -fc_threshold), 'significance'] = 'Downregulated'

    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=de_results, 
        x='log2FC', 
        y=-np.log10(de_results['adj_p_value']),
        hue='significance',
        palette={'Upregulated': 'red', 'Downregulated': 'blue', 'Not Significant': 'grey'},
        alpha=0.6,
        s=20
    )
    
    # Add annotation for top genes
    top_hits = de_results[de_results['significance'] != 'Not Significant'].sort_values('adj_p_value').head(10)
    for i, row in top_hits.iterrows():
        plt.text(row['log2FC'], -np.log10(row['adj_p_value']), i, fontsize=9)

    plt.title('Volcano Plot', fontsize=16)
    plt.xlabel('Log2 Fold Change', fontsize=12)
    plt.ylabel('-log10(Adjusted P-value)', fontsize=12)
    plt.axvline(fc_threshold, color='black', linestyle='--', lw=1)
    plt.axvline(-fc_threshold, color='black', linestyle='--', lw=1)
    plt.axhline(-np.log10(p_threshold), color='black', linestyle='--', lw=1)
    plt.legend(title='Significance')
    plt.grid(alpha=0.3)
    # Save the figure to be used in the README
    plt.savefig('../assets/volcano_plot_example.png', dpi=300, bbox_inches='tight')
    plt.show()
