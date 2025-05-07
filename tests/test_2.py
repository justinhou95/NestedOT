def test_threads():
    # Test the number of threads does not affect
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

    markovian = True
    num_threads_list = [0, 1, 8, 16]
    for markovian in [True, False]:
        for p in [0, 1, 2, 4]:
            print("Markovian: ", markovian, " p: ", p)
            v_list = []
            for num_threads in num_threads_list:
                print("NumThread: ", num_threads)
                v = nested_ot(
                    X,
                    Y,
                    grid_size,
                    markovian,
                    num_threads=num_threads,
                    power=2,
                    verbose=True,
                )
                v_list.append(v)
            for v in v_list:
                assert np.abs(v_list[0] - v) < 1e-10
            print("=" * 100)


if __name__ == "__main__":
    test_threads()
