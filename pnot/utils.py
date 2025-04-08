import matplotlib.pyplot as plt
import numpy as np


def Lmatrix2paths(L, n_sample, normalize=False, seed=0, verbose=False):
    r"""
    Lower triangular matrix L to covariance matrix A and generated paths
    """
    A0 = L @ L.T  # A = LL^T
    L = L / np.sqrt(np.trace(A0)) if normalize else L
    A = L @ L.T

    # np.linalg.cholesky(A) -  L (sanity check)

    if verbose:
        print("Cholesky:")
        print(L)
        print("Covariance:")
        print(A)

    T = len(L)

    np.random.seed(seed)
    noise1 = np.random.normal(size=[T, n_sample])  # (T, n_sample)
    X = L @ noise1  # (T, n_sample)
    X = np.concatenate([np.zeros_like(X[:1]), X], axis=0)  # (T+1, n_sample)
    return X, A


def nested(X, Y):
    import _wrapper

    return _wrapper.nested(X, Y)
