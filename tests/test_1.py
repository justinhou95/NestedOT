def test_1():
    import numpy as np
    from pnot import nested_ot
    from pnot.utils import adapted_wasserstein_squared, matrixL2paths

    random_seed = 0
    n_sample = 100
    T = 3

    L = np.array([[1, 0, 0], [1, 1, 0], [1, 1, 1]])
    X, A = matrixL2paths(L, n_sample, seed=random_seed)

    M = np.array([[1, 0, 0], [2, 1, 0], [2, 1, 2]])
    Y, B = matrixL2paths(M, n_sample, seed=random_seed)

    awd2square = adapted_wasserstein_squared(A, B)

    grid_size = 0.1

    from pnot.py_solver import nested_ot_solver_py

    for markovian in [True, False]:
        v1 = nested_ot_solver_py(X, Y, grid_size, markovian, parallel=True)
        v2 = nested_ot_solver_py(X, Y, grid_size, markovian, parallel=True)
        v3 = nested_ot(X, Y, grid_size, markovian)
        print(np.abs(v1 - v2))
        print(np.abs(v2 - v3))
        assert np.abs(v1 - v2) < 1e-10
        assert np.abs(v2 - v3) < 1e-10


if __name__ == "__main__":
    test_1()
