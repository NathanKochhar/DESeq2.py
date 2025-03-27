'''
We want to normalize the counts data due to sequencing depth. Samples might have a different
number of total counts which may cause downstream problems. We will normalize to scale the counts
so they are even across the board.

Steps are:
- Calculating the mean of counts for each gene
- Dividing the counts of each gene by its geometric mean across samples
- Taking the median of these ratios to derive size factors that are used to scale the counts.
'''

import pandas as pd
import numpy as np

def normalization(raw_counts):

    # The geometric mean of these counts is computed for each gene. 
    # This is achieved using the natural logarithm to mitigate extreme values and then exponentiating the average of these logarithms. 
    # This step helps standardize counts for genes across samples.
    geometric_means = np.exp(np.log(raw_counts).mean(axis=1))
    print(geometric_means)

    # Each gene's count in the data frame is divided by its geometric mean.
    # This division adjusts the counts according to the central tendency specific to the gene.
    size_factors = raw_counts.divide(geometric_means, axis=0)
    print(size_factors)

    return size_factors.median(axis=1)
