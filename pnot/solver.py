from collections import defaultdict
import concurrent.futures
import copy
import random

import concurrent
import matplotlib.pyplot as plt
import numpy as np
import ot
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm


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
    """
    Project paths to adapted grids
    """
    grid_func = lambda x: np.floor(x / delta_n + 0.5) * delta_n
    adapted_samples = grid_func(samples)
    return adapted_samples


def sort_qpath(path):
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


class ConditionalLaw:
    r"""For easy notation and structure"""

    def __init__(self, qX, markovian=True):
        self.T = qX.shape[-1] - 1
        self.markovian = markovian
        self.mu_x = qpath2mu_x(qX, markovian)
        self.c, self.nc, self.v, self.w, self.next_idx, self.nv_cum, self.q2idx = (
            list_repr_mu_x(self.mu_x, markovian)
        )


def solve_ot(wx, wy, cost):
    if len(wx) == 1 or len(wy) == 1:
        res = np.dot(np.dot(wx, cost), wy)  # in this case we has closed solution
    else:
        res = np.sum(
            cost * ot.lp.emd(wx, wy, cost)
        )  # faster than ot.emd2(wx, wy, cost)
    return res


def nested2(kernel_x: ConditionalLaw, kernel_y: ConditionalLaw, cost_matrix):
    assert kernel_x.T == kernel_y.T
    T = kernel_x.T
    V = [
        np.zeros([kernel_x.nc[t], kernel_y.nc[t]]) for t in range(T)
    ]  # V_t(x_{1:t},y_{1:t})
    for t in range(T - 1, -1, -1):
        x_bar = tqdm(range(kernel_x.nc[t]))
        x_bar.set_description(f"Timestep {t}")
        for ix in x_bar:
            for iy in range(kernel_y.nc[t]):
                vx = kernel_x.v[t][ix]
                vy = kernel_y.v[t][iy]
                wx = kernel_x.w[t][ix]
                wy = kernel_y.w[t][iy]
                cost = cost_matrix[np.ix_(vx, vy)]
                if t < T - 1:
                    x_next_idx = kernel_x.next_idx[t][ix]
                    y_next_idx = kernel_y.next_idx[t][iy]
                    cost += V[t + 1][np.ix_(x_next_idx, y_next_idx)]
                V[t][ix, iy] = solve_ot(wx, wy, cost)
    AW_2square = V[0][0, 0]
    return AW_2square


def chunk_process(arg):
    x_arg, y_arg, Vtplus, cost_matrix = arg
    x_arg[0] = tqdm(x_arg[0])
    Vt = np.zeros([len(x_arg[0]), len(y_arg[0])])
    if Vtplus is None:
        for cx, vx, wx in zip(*x_arg[:-1]):
            for cy, vy, wy in zip(*y_arg[:-1]):
                cost = cost_matrix[np.ix_(vx, vy)]
                Vt[cx, cy] = solve_ot(wx, wy, cost)
    else:
        for cx, vx, wx, idx_x in zip(*x_arg):
            for cy, vy, wy, idx_y in zip(*y_arg):
                cost = cost_matrix[np.ix_(vx, vy)]
                cost += Vtplus[np.ix_(idx_x, idx_y)]
                Vt[cx, cy] = solve_ot(wx, wy, cost)
    return Vt


def nested2_parallel(kernel_x: ConditionalLaw, kernel_y: ConditionalLaw, cost_matrix):
    r"""This only works wtih non-Markovian now because slicing Vtplus is easier with non-Markovian"""
    assert kernel_x.T == kernel_y.T

    T = kernel_x.T
    V = [
        np.zeros([kernel_x.nc[t], kernel_y.nc[t]]) for t in range(T)
    ]  # V_t(x_{1:t},y_{1:t})
    for t in range(T - 1, -1, -1):
        n_processes = 8 if t > 0 else 1  # HERE WE NEED TO CHANGE BACK TO t>1
        chunks = np.array_split(range(kernel_x.nc[t]), n_processes)
        args = []
        for chunk in chunks:
            ix_start, ix_end = chunk[0], chunk[-1] + 1
            x_arg = [
                range(len(chunk)),
                kernel_x.v[t][ix_start:ix_end],
                kernel_x.w[t][ix_start:ix_end],
                kernel_x.next_idx[t][ix_start:ix_end],
            ]
            y_arg = [
                range(kernel_y.nc[t]),
                kernel_y.v[t],
                kernel_y.w[t],
                kernel_y.next_idx[t],
            ]
            Vtplus = V[t + 1] if t < T - 1 else None
            arg = (x_arg, y_arg, Vtplus, cost_matrix)
            args.append(arg)

        with concurrent.futures.ProcessPoolExecutor() as executor:
            Vts = executor.map(chunk_process, args)
        for chunk, Vt in zip(chunks, Vts):
            V[t][chunk] = Vt

        # for arg, chunk in zip(args, chunks):
        #     res = chunk_process(arg)
        #     V[t][chunk] = res

    AW_2square = V[0][0, 0]
    return AW_2square
