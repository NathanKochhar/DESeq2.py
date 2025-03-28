'''
We want to normalize the counts data due to sequencing depth. Samples might have a different
number of total counts which may cause downstream problems. We will normalize to scale the counts
so they are even across the board.

Steps are:
- Calculating the geometric mean of counts for each gene
- Dividing the counts of each gene by its geometric mean across samples
- Taking the median of these ratios by sample to derive size factors
- Divide raw counts by sample size factors to get the normalized matrix
'''

import pandas as pd
import numpy as np

def calculate_geometric_means(raw_counts: pd.DataFrame):
    log_df = np.log(raw_counts.replace(0, np.nan))
    geo_means = np.exp(log_df.mean(axis=1, skipna=True))
    return(geo_means)

def calculate_ratios(raw_counts: pd.DataFrame, geo_means: pd.Series):
    ratios = raw_counts.divide(geo_means, axis=0)
    return(ratios)

def estimate_size_factors(ratios: pd.DataFrame):
    size_factors = ratios.median(axis=0)
    return(size_factors)

def normalize_counts(raw_counts: pd.DataFrame, size_factors: pd.Series):
    norm_counts = raw_counts.divide(size_factors, axis=1)
    return(norm_counts)