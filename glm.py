'''
glm.py

'''

import pandas as pd
import numpy as np
import math

def build_design_matrix(metadata: pd.DataFrame, formula: str = "~ condition"):
    variable_name = formula.split("~")[1].strip()
    print(metadata[variable_name])
    metadata = pd.get_dummies(metadata, columns=[variable_name], drop_first=True, dtype='int')
    metadata['intercept'] = 1
    metadata2 = metadata.iloc[:, [-1,1]]
    return(metadata2)

def initialize_mu(norm_counts: pd.DataFrame, size_factors: pd.Series):
    means = norm_counts.mean(axis=1)
    mu_mtx = pd.DataFrame(index=norm_counts.index, columns=norm_counts.columns)
    for i in range(0, int(mu_mtx.shape[0])): #row
        for j in range(0, (int(mu_mtx.shape[1]))): #column
            mu_mtx.iloc[i,j] = means.iloc[i] * size_factors.iloc[j]
            #print("i: ",i,", j: ", j,", mean: ", means.iloc[i],", size_factors: ", size_factors.iloc[j], "mu: ", means.iloc[i] * size_factors.iloc[j])
    return(mu_mtx)


def calculate_weights(mu: pd.DataFrame, dispersions: pd.Series):
    alpha = dispersions
    alpha = alpha.reindex(mu.index)
    mu_values = mu.values
    alpha_values = alpha.values[:, np.newaxis]
    variance = mu_values + alpha_values * mu_values**2
    weights = mu_values**2 / variance
    return(pd.DataFrame(weights, index=mu.index, columns=mu.columns))

def check_convergence(old_mu: pd.DataFrame, new_mu: pd.DataFrame, tol: float):
    eps = 1e-8
    old_mu = old_mu[new_mu.columns]
    old_mu = old_mu.loc[new_mu.index]
    relative_change = (new_mu - old_mu).abs() / (old_mu.abs() + eps)
    max_change = relative_change.values.max()
    return max_change < tol


def fit_glm_nb(counts: pd.DataFrame, design_matrix: pd.DataFrame, dispersions: pd.Series, size_factors: pd.Series, max_iter: int=50, tol: float=1e-6):
    initial_mu = initialize_mu(counts, size_factors)
    print(initial_mu)
    linear_predictor = np.log10(initial_mu.astype(np.float64))
    print(linear_predictor)
    beta = pd.DataFrame(data = 0.0, index=counts.index, columns=design_matrix.columns)
    print(beta)

    '''genes = counts.index
    samples = counts.columns
    predictors = design_matrix.columns

    betas = pd.DataFrame(index=genes, columns=predictors, dtype=float)
    mu = pd.DataFrame(index=genes, columns=samples, dtype=float)
    convergence = pd.Series(False, index=genes)

    for gene in genes:
        y = counts.loc[gene]
        alpha = dispersions[gene]
        s = size_factors

        # Initial mu and eta
        y_norm = y / s
        mu_g = s * y_norm.mean()
        eta = np.log(mu_g + 1e-8)  # avoid log(0)
        beta = np.zeros(len(predictors))

        converged = False

        for _ in range(max_iter):
            # 1. Compute weights
            mu_df = pd.DataFrame([mu_g], index=[gene], columns=samples)
            weights = calculate_weights(mu_df, pd.Series({gene: alpha})).loc[gene]

            # 2. Compute working response z
            dmu_deta = mu_g
            z = eta + (y - mu_g) / (dmu_deta + 1e-8)

            # 3. Weighted least squares
            W = np.diag(weights)
            X = design_matrix.values
            XtW = X.T @ W
            try:
                beta_new = np.linalg.solve(XtW @ X, XtW @ z)
            except np.linalg.LinAlgError:
                break  # Singular matrix — skip this gene

            # 4. Update eta and mu
            eta_new = X @ beta_new
            mu_new = s * np.exp(eta_new)

            # 5. Check convergence
            mu_old_df = pd.DataFrame([mu_g], index=[gene], columns=samples)
            mu_new_df = pd.DataFrame([mu_new], index=[gene], columns=samples)
            if check_convergence(mu_old_df, mu_new_df, tol):
                converged = True
                break

            # 6. Update for next iteration
            beta = beta_new
            mu_g = mu_new
            eta = eta_new

        # Store results
        betas.loc[gene] = beta
        mu.loc[gene] = mu_g
        convergence[gene] = converged

    return betas, mu, convergence'''





    return("hi")
    #should return : Tuple[betas, mu, converged]