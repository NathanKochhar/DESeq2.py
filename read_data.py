'''
Functions to read in the counts matrix and metadata
'''

import pandas as pd

def read_mtx(path):
    if path.split(".")[-1] == "csv":
        return(pd.read_csv(path, index_col=0))
    if path.split(".")[-1] == "tsv":
        return(pd.read_tsv(path, index_col=0, sep='\t'))
    
def read_conds(path):
    return(pd.read_csv(path, index_col=0))






