'''
estimate dispersion
'''

import pandas as pd
import numpy as np

def calculate_means(norm_counts: pd.DataFrame):
    means = norm_counts.mean(axis=1)
    return(means)

def calculate_variances(norm_counts: pd.DataFrame):
    variances = norm_counts.var(axis=1)
    return(variances)

def estimate_raw_dispersion(means: pd.Series, variances: pd.Series):
    raw_dispersion = (variances*means)/(means*means)
    return(raw_dispersion)

def filter_genes_for_dispersion(norm_counts: pd.DataFrame, min_counts: int = 10):
    norm_counts['row_sum'] = norm_counts.sum(axis=1)
    filtered_df = norm_counts[norm_counts['row_sum'] > min_counts]
    return(filtered_df)

def get_dispersion_table(raw_dispersion: pd.Series, means: pd.Series):
    data = {'raw_dispersion': raw_dispersion, 'means': means}
    df = pd.DataFrame(data)
    return(df)

