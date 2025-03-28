'''
Functions to read in the counts matrix and metadata
'''

import pandas as pd

def load_counts(path: str):
    if path.split(".")[-1] == "csv":
        return(pd.read_csv(path, index_col=0))
    if path.split(".")[-1] == "tsv":
        return(pd.read_tsv(path, index_col=0, sep='\t'))
    
def load_metadata(path: str):
    return(pd.read_csv(path, index_col=0))






