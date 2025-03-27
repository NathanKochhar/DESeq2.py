'''
This is where we will run the workflow.

'''

import read_data as rd
import normalization as nm

Counts_Matrix_Path = "example_data/counts_matrix.csv"
data = rd.read_mtx(Counts_Matrix_Path)

#remove rows with rowsum less than 10
row_sums = data.sum(axis=1)
rows_to_remove = row_sums[row_sums < 10].index
data = data.drop(rows_to_remove)
print(data.head())

Conds_Path = "example_data/conds.csv"
conds = rd.read_conds(Conds_Path)
#print(conds.head())

normalized = nm.normalization(data)

print(normalized)
