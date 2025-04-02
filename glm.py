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

'''
/Users/nathanielkochhar/Documents/DESeq2.py/glm.py:33: RuntimeWarning: overflow encountered in square
  variance = mu_values + alpha_values * mu_values**2
/Users/nathanielkochhar/Documents/DESeq2.py/glm.py:34: RuntimeWarning: overflow encountered in square
  weights = mu_values**2 / variance
/Users/nathanielkochhar/Documents/DESeq2.py/glm.py:34: RuntimeWarning: invalid value encountered in divide
  weights = mu_values**2 / variance
'''

def check_convergence(old_mu: pd.DataFrame, new_mu: pd.DataFrame, tol: float):
    eps = 1e-8
    old_mu = old_mu[new_mu.columns]
    old_mu = old_mu.loc[new_mu.index]
    relative_change = (new_mu - old_mu).abs() / (old_mu.abs() + eps)
    max_change = relative_change.values.max()
    return max_change < tol


def fit_glm_nb(counts: pd.DataFrame, design_matrix: pd.DataFrame, dispersions: pd.Series, size_factors: pd.Series, max_iter: int=50, tol: float=1e-6):
    
    '''
    initial_mu = initialize_mu(counts, size_factors)
    print(initial_mu)
    linear_predictor = np.log10(initial_mu.astype(np.float64))
    print(linear_predictor)
    beta = pd.DataFrame(data = 0.0, index=counts.index, columns=design_matrix.columns)
    print(beta)
    '''


    genes = counts.index
    samples = counts.columns
    predictors = design_matrix.columns

    betas = pd.DataFrame(data = 0.0, index=counts.index, columns=design_matrix.columns)
    mu = initialize_mu(counts, size_factors)
    convergence = pd.Series(False, index=genes)

    for gene in genes:
        #print("\n",gene, "\n")
        y = counts.loc[gene]
        alpha = dispersions[gene]
        s = size_factors

        # Initial mu and eta
        '''y_norm = y / s
        mu_g = s * y_norm.mean()'''
        mu_g = mu.loc[gene]
        eta = np.log10(mu_g.astype(np.float64)+ 1e-8)  # avoid log(0)     linear_predictor
        #np.log10(mu_g.astype(np.float64)+ 1e-8)
        beta = np.zeros(len(predictors))

        converged = False
        #print("\n gene: ",gene, "\n")

        for _ in range(max_iter):
            # 1. Compute weights
            mu_df = pd.DataFrame([mu_g], index=[gene], columns=samples)
            weights = calculate_weights(mu_df, pd.Series({gene: alpha})).loc[gene]
            #print("weights: ",weights)

            # 2. Compute working response z
            dmu_deta = mu_g
            z = eta + (y - mu_g) / (dmu_deta + 1e-8)
            #print("\n working response z: ",z)

            # 3. Weighted least squares
            W = np.diag(weights)
            #print("\n W: ",W)
            X = design_matrix.values
            #print("\n X.T: ",X.T)
            XtW = X.T @ W
            #print("\n XtW: ",XtW)

            try:
                #beta_new = np.linalg.solve(XtW @ X, XtW @ z)
                XtWX = (XtW @ X).astype(np.float64)
                XtWz = (XtW @ z).astype(np.float64)
                beta_new = np.linalg.solve(XtWX, XtWz)

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
                print("gene: ",gene, "- Converged on iteration: ", _ + 1)
                break

            #print("gene: ",gene, "- iteration: ", _ + 1)

            # 6. Update for next iteration
            beta = beta_new
            mu_g = mu_new
            eta = eta_new
        if converged == False:
            print("gene: ",gene, "- didn't converge")
        # Store results
        betas.loc[gene] = beta
        mu.loc[gene] = mu_g
        convergence[gene] = converged

    return betas, mu, convergence


def calculate_standard_errors(design_matrix: pd.DataFrame, weights: pd.DataFrame, mu: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate standard errors for NB GLM coefficients using the Fisher Information matrix.

    Parameters:
    -----------
    design_matrix : pd.DataFrame
        Design matrix (samples x predictors).
    weights : pd.DataFrame
        Observation-level weights (genes x samples).
    mu : pd.DataFrame
        Fitted means (genes x samples).

    Returns:
    --------
    pd.DataFrame
        Standard errors for each coefficient (genes x predictors).
    """
    genes = weights.index
    predictors = design_matrix.columns
    samples = design_matrix.index

    # Output SEs DataFrame
    se = pd.DataFrame(index=genes, columns=predictors, dtype=float)

    X = design_matrix.loc[samples].values.astype(np.float64)

    for gene in genes:
        # Get weights for this gene
        w = weights.loc[gene].reindex(samples).values.astype(np.float64)
        
        # Form diagonal weight matrix W_i
        W = np.diag(w)
        
        try:
            # Compute Fisher Information Matrix: XtWX = X^T W X
            XtWX = X.T @ W @ X

            # Invert to get covariance matrix
            XtWX_inv = np.linalg.inv(XtWX)

            # Standard errors = sqrt of diagonal of covariance matrix
            se_values = np.sqrt(np.diag(XtWX_inv))

            # Store in output DataFrame
            se.loc[gene] = se_values

        except np.linalg.LinAlgError:
            # If singular matrix, fill with NaNs for that gene
            se.loc[gene] = np.nan

    return se
