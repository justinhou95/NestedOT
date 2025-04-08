import numpy as np
from pnot.utils import Lmatrix2paths, nested

if __name__ == "__main__":

    random_seed = 0
    n_sample = 1000
    T = 3
    L = np.array([[1, 0, 0], [2, 4, 0], [3, 2, 1]])
    X, A = Lmatrix2paths(L, n_sample, seed=random_seed, verbose=False)
    M = np.array([[1, 0, 0], [2, 3, 0], [3, 1, 2]])
    Y, B = Lmatrix2paths(M, n_sample, seed=random_seed, verbose=False)
    print(L)
    print(M)

    nested(X, Y)
