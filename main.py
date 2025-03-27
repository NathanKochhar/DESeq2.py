'''
This is where we will run the workflow.

'''

import read_data as rd

Counts_Matrix_Path = "example_data/counts_matrix.csv"
data = rd.read_mtx(Counts_Matrix_Path)
print(data.head())

Conds_Path = "example_data/conds.csv"
conds = rd.read_conds(Conds_Path)
print(conds.head())


