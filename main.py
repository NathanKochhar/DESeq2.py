'''
This is where we will run the workflow.

'''

import read_data as rd
import normalization as nm
import estimate_dispersion as ed

Counts_Matrix_Path = "example_data/counts_matrix.csv"
data = rd.load_counts(Counts_Matrix_Path)
print("THIS IS THE RAW MATRIX: \n")
#print(data)

Conds_Path = "example_data/conds.csv"
conds = rd.load_metadata(Conds_Path)
#print(conds.head())

gm = nm.calculate_geometric_means(data)
print("THESE ARE MY GEOMETRIC MEANS FOR EACH GENE \n")
#print(gm)

ratios = nm.calculate_ratios(data, gm)
print("THESE ARE MY RATIOS \n")
#print(ratios)

size_factors = nm.estimate_size_factors(ratios)
print("THESE ARE MY SIZE FACTORS \n")
#print(size_factors)

norm_counts = nm.normalize_counts(data, size_factors)
print("THIS IS MY NORMALIZED MATRIX \n")
print(norm_counts)

means = ed.calculate_means(norm_counts)
print("THESE ARE MY GENE MEANS \n")
#print(means)

variances = ed.calculate_variances(norm_counts)
print("THESE ARE MY GENE VARIANCES \n")
#print(variances)

raw_dispersion = ed.estimate_raw_dispersion(means, variances)
print("THESE ARE MY RAW DISPERSIONS \n")
#print(raw_dispersion)

filtered_matrix = ed.filter_genes_for_dispersion(norm_counts, 5)
print("THIS IS MY FILTERED NORMALIZED MATRIX \n")
print(filtered_matrix)

dispersion_table = ed.get_dispersion_table(raw_dispersion, means)
print("THIS IS MY DISTPERSION TABLE \n")
print(dispersion_table)
