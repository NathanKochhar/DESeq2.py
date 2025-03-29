'''
negative_binomial.py

'''

import pandas as pd
import numpy as np
import math

def nb_variance(mu: np.ndarray, alpha: np.ndarray):
    nb_variance = (mu) + (alpha*(mu**2))
    return (nb_variance)

#Stringlings approximation
def log_gamma(x: float):
    if (x <= 0):
        raise Exception("Sorry, you can't log_gamma a 0")
    elif (x < 1):
        x = x + 1
        log_gamma = x * (np.log10(x)) - x + (.5*np.log10(2*np.pi)) + (.5*np.log10(x)) - np.log10(x-1)
        return(log_gamma)
    else:
        log_gamma = x * (np.log10(x)) - x + (.5*np.log10(2*np.pi)) + (.5*np.log10(x))
        return(log_gamma)

'''
shape[1] is cols
Can do for i in range(0:(matrix.shape[1]-1))

'''
def nb_log_pdf(matrix: np.ndarray, mu: np.ndarray, alpha: np.ndarray):
    log_likelihood_mtx = pd.DataFrame(index=matrix.index, columns=matrix.columns)
    num_cols = int(matrix.shape[1]) 
    num_rows = int(matrix.shape[0])
    for i in range(0, num_rows): #row?
        for j in range(0, num_cols): #column?
            y = matrix.iloc[i,j]
            mean = mu.iloc[i]
            disp = alpha.iloc[i]
            #print("row: ",i+1,", column: ", j+1,", value: ", y,", mean: ", mean,", disp: ", disp)
            if y < 20:
                log_fact = np.log10((math.factorial(int(y))))
            else:
                log_fact = 0.5 * np.log10(2*np.pi*y) + y*np.log10(y/np.e)
            nb_log_likelihood_val_1 = log_gamma(y+(1/disp)) - log_gamma(1/disp) - log_fact + (1/disp)*np.log10(1/(1+(mean*disp)))
            nb_log_likelihood_val = nb_log_likelihood_val_1 + y*np.log10((mean*disp)/(1+(mean*disp)))
            log_likelihood_mtx.iloc[i, j] = nb_log_likelihood_val
    #print(log_likelihood_mtx)
    return(log_likelihood_mtx)
