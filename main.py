'''
This is where we will run the workflow.

'''

import read_data as rd
import normalization as nm
import estimate_dispersion as ed
import negative_binomial as nb
import glm

Counts_Matrix_Path = "example_data/counts_matrix.csv"
data = rd.load_counts(Counts_Matrix_Path)
data = data.iloc[:, 5:11]
print("\n THIS IS THE RAW MATRIX: \n")
#print(data)

Conds_Path = "example_data/conds.csv"
conds = rd.load_metadata(Conds_Path)
conds = conds.iloc[5:11, :]
print("\n THIS IS THE CONDS TABLE: \n")
print(conds)

norm_counts, size_factors = nm.run_normalization(data)
print("\n THIS IS MY NORMALIZED MATRIX: \n")
#print(norm_counts)

filtered_matrix, dispersion_table = ed.run_dispersion(norm_counts, min_counts = 5)
print("\n THIS IS MY FILTERED NORMALIZED MATRIX: \n")
print(filtered_matrix)
print("\n THIS IS MY DISTPERSION TABLE: \n")
print(dispersion_table)

means = dispersion_table['means']
raw_dispersion = dispersion_table['raw_dispersion']
nb_var = nb.nb_variance(means, raw_dispersion)
print("\n THIS IS MY NEGITIVE BIONOMIAL VARIANCE: \n")
#print(nb_var)

'''nb_log_pdf = nb.nb_log_pdf(filtered_matrix, means, raw_dispersion)
print("\n THIS IS MY NEGITIVE BIONOMIAL LOG LIKELYHOOD PDF: \n")
print(nb_log_pdf)'''

design_matrix = glm.build_design_matrix(conds)
print("\n THIS IS MY DESIGN MATRIX: \n")
print(design_matrix)

'''initialize_mu = glm.initialize_mu(filtered_matrix, size_factors)
print("\n THIS IS MY INITIALIZED MU: \n")
print(initialize_mu)'''

'''print("\n CALC WEIGHTS TEST: \n")
weights = glm.calculate_weights(initialize_mu, raw_dispersion)
print(weights)'''

fit_glm = glm.fit_glm_nb(filtered_matrix, design_matrix, raw_dispersion, size_factors, max_iter=50, tol=1e-5)
print("\n THIS IS MY FIT NB GLM: \n")
print(fit_glm)



