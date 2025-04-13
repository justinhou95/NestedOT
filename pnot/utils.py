from collections import defaultdict

import numpy as np


def matrixL2paths(L, n_sample, seed=0):
    r"""
    Lower triangular matrix L to covariance matrix A and generated paths
    """
    A = L @ L.T
    T = len(L)
    np.random.seed(seed)
    noise1 = np.random.normal(size=[T, n_sample])  # (T, n_sample)
    X = L @ noise1  # (T, n_sample)
    X = np.concatenate([np.zeros_like(X[:1]), X], axis=0)  # (T+1, n_sample)
    return X, A


def adapted_wasserstein_squared(A, B, a=0, b=0):
    # Cholesky decompositions: A = L L^T, B = M M^T
    L = np.linalg.cholesky(A)
    M = np.linalg.cholesky(B)
    # Mean squared difference
    mean_diff = np.sum((a - b) ** 2)
    # Trace terms
    trace_sum = np.trace(A) + np.trace(B)
    # L1 norm of diagonal elements of L^T M
    l1_diag = np.sum(np.abs(np.diag(L.T @ M)))
    # Final adapted Wasserstein squared distance
    return mean_diff + trace_sum - 2 * l1_diag


def path2adaptedpath(samples, delta_n):
    r"""
    Project paths to adapted grids
    """
    grid_func = lambda x: np.floor(x / delta_n + 0.5) * delta_n
    adapted_samples = grid_func(samples)
    return adapted_samples


def sort_qpath(path):
    r"""
    Sort the quantized paths iteratively by xq_0, xq_1, xq_2, ...
    """
    T = path.shape[-1] - 1
    sorting_keys = [path[:, i] for i in range(T, -1, -1)]
    return path[np.lexsort(tuple(sorting_keys))]


def qpath2mu_x(qpath, markovian=False):
    r"""
    Quantized Path to Conditional Measure
    non-Markovian:
    mu_x[0] = {(3,): {1: 1, 2: 5}}
    Markovian:
    mu_x[0] = {3: {1: 1, 2: 5}}
    """
    T = qpath.shape[-1] - 1
    mu_x = [defaultdict(dict) for t in range(T)]
    for t in range(T):
        for path in qpath:
            if markovian:
                pre_path = (path[t],)
            else:
                pre_path = tuple([x for x in path[: t + 1]])
            next_val = int(path[t + 1])
            if pre_path not in mu_x[t] or next_val not in mu_x[t][pre_path]:
                mu_x[t][pre_path][next_val] = 1
            else:
                mu_x[t][pre_path][next_val] += 1
    return mu_x


def list_repr_mu_x(mu_x, markovian=False):
    r"""
    represent mu_x[t] with
    mu_x_c[t][i]: xq_{1:t} quantized conditional path up to time t
    mu_x_v[t][i]: a list of values xq_{t+1} follows xq_{1:t}
    mu_x_w[t][i]: a list of weights mu_{x_{1:t}}(x_{t+1})

    e.g. if we have 4 paths (same paths count twice)
    quantized
    (1, 2, 3)
    (1, 2, 4)
    (1, 2, 4)
    (2, 3, 5)
    Then we have:
    mu_x_c[t=1] = [(1,2), (2,3)]
    mu_x_v[t=1][i=1] = [3, 4]
    mu_x_w[t=1][i=1] = [1/3, 2/3]

    """
    T = len(mu_x)
    mu_x_c = [None for _ in range(T)]  # mu_x_c[t][ix] = xq_{1:t}
    mu_x_nc = [None for _ in range(T)]  # mu_x_nc[t] = number of xq_{1:t}
    mu_x_v = [[] for _ in range(T)]  # mu_x_v[t][ix] = a list of values of xq_{t+1}
    mu_x_w = [[] for _ in range(T)]  # mu_x_w[t][ix] = a list of weights of xq_{t+1}
    for t in range(T - 1, -1, -1):
        mu_x_t = mu_x[t]
        # conditions
        pre_paths = list(mu_x_t.keys())
        mu_x_c[t] = pre_paths
        mu_x_nc[t] = len(pre_paths)
        # distributions
        for dist in mu_x_t.values():
            # values
            values = np.array(list(dist.keys()))
            mu_x_v[t].append(values)
            # weights
            counts = np.array(list(dist.values()))
            weights = counts / np.sum(counts)
            mu_x_w[t].append(weights)

    mu_x_nv_cum = [
        [] for _ in range(T)
    ]  # mu_x_nv_cum[t] = cumsum list of len(mu_x_v[t][ix])
    mu_x_q2idx = [
        None for _ in range(T)
    ]  # mu_x_q2idx[t] = map from mu_x_c[t][ix] to ix
    mu_x_next_idx = [
        [] for _ in range(T)
    ]  # mu_x_next_idx[t][ix] = a list of index in Vtplus of xq_{t+1}

    for t in range(T - 1, -1, -1):
        mu_x_nv = [len(vs) for vs in mu_x_v[t]]
        mu_x_nv_cum[t] = np.cumsum([0] + mu_x_nv)
        mu_x_q2idx[t] = {c[-1]: i for i, c in enumerate(mu_x_c[t])}
        if t < T - 1:
            for ix in range(mu_x_nc[t]):
                if markovian:
                    mu_x_next_idx[t].append(
                        [mu_x_q2idx[t + 1][v] for v in mu_x_v[t][ix]]
                    )
                else:
                    mu_x_next_idx[t].append(
                        range(mu_x_nv_cum[t][ix], mu_x_nv_cum[t][ix + 1])
                    )

    return mu_x_c, mu_x_nc, mu_x_v, mu_x_w, mu_x_next_idx, mu_x_nv_cum, mu_x_q2idx
