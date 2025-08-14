import pandas as pd
import numpy as np

def load_and_clean_maxquant(protein_groups_path, metadata_path):
    """
    Loads a MaxQuant proteinGroups.txt file, filters it, and normalizes intensities.
    """
    # Load data
    df = pd.read_csv(protein_groups_path, sep='\t')
    meta = pd.read_csv(metadata_path)

    # 1. Filter out contaminants and reverse hits
    df = df[df['Reverse'] != '+']
    df = df[df['Potential contaminant'] != '+']

    # 2. Extract intensity columns and gene names
    intensity_cols = [f'Intensity {s}' for s in meta['sample']]
    df_intensity = df[['Gene names'] + intensity_cols].copy()
    df_intensity = df_intensity.dropna(subset=['Gene names'])
    df_intensity = df_intensity.drop_duplicates(subset=['Gene names'], keep='first')
    df_intensity = df_intensity.set_index('Gene names')

    # 3. Replace 0 with NaN and log2 transform
    df_intensity[df_intensity == 0] = np.nan
    df_log2 = np.log2(df_intensity)

    # 4. Median normalization
    median_diff = df_log2.median(axis=0) - df_log2.median(axis=0).median()
    df_normalized = df_log2 - median_diff

    # 5. Impute missing values (simple mean imputation)
    df_imputed = df_normalized.apply(lambda x: x.fillna(x.mean()), axis=1)

    return df_imputed, df_imputed.index.tolist()
