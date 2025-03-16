# Define your function for Task 2 here.
import numpy as np

def fd_coeffs(beta):
    '''
    Compute finite difference coefficients for a given set of beta values.
    Parameters
    ----------
    beta : array_like
        Array of beta values.
    Returns
    -------
    alpha : ndarray
        Array of finite difference coefficients.
    '''
    N = len(beta)
    # initialize the matrix and the coefficient array
    A = np.ones((N, N))
    alpha = np.zeros(N)
    y = np.zeros(N)
    # let the coefficient of the first derivative be 1
    y[1] = 1
    # add the coefficients of the each order of derivatives to the matrix
    for i in range(1,N):
        A[i, :] = beta ** i/ np.math.factorial(i)
    # use linear solving to find the coefficients  
    alpha = np.linalg.solve(A, y)
    return alpha