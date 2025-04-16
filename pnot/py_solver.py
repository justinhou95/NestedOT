import concurrent
import concurrent.futures
import time

import numpy as np
import ot
from tqdm import tqdm

from .utils import list_repr_mu_x, path2adaptedpath, qpath2mu_x, sort_qpath


class ConditionalLaw:
    r"""For easy notation and tructure"""

    def __init__(self, qX, markovian=True):
        self.T = qX.shape[-1] - 1
        self.markovian = markovian
        self.mu_x = qpath2mu_x(qX, markovian)
        self.c, self.nc, self.v, self.w, self.next_idx, self.nv_cum, self.q2idx = (
            list_repr_mu_x(self.mu_x, markovian)
        )


def ot_solver(wx, wy, cost):
    if len(wx) == 1 or len(wy) == 1:
        opt_value = np.dot(np.dot(wx, cost), wy)  # in this case we has closed solution
    else:
        opt_cpl = ot.lp.emd(wx, wy, cost)
        opt_value = np.sum(cost * opt_cpl)  # faster than ot.emd2(wx, wy, cost)
    return opt_value


def nested(kernel_x: ConditionalLaw, kernel_y: ConditionalLaw, cost_matrix):
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
                V[t][ix, iy] = ot_solver(wx, wy, cost)
    nested_ot_value = V[0][0, 0]
    return nested_ot_value


def chunk_process(arg):
    x_arg, y_arg, Vtplus, cost_matrix = arg
    x_arg[0] = tqdm(x_arg[0])
    Vt = np.zeros([len(x_arg[0]), len(y_arg[0])])
    if Vtplus is None:
        for cx, vx, wx in zip(*x_arg[:-1]):
            for cy, vy, wy in zip(*y_arg[:-1]):
                cost = cost_matrix[np.ix_(vx, vy)]
                Vt[cx, cy] = ot_solver(wx, wy, cost)
    else:
        for cx, vx, wx, idx_x in zip(*x_arg):
            for cy, vy, wy, idx_y in zip(*y_arg):
                cost = cost_matrix[np.ix_(vx, vy)]
                cost += Vtplus[np.ix_(idx_x, idx_y)]
                Vt[cx, cy] = ot_solver(wx, wy, cost)
    return Vt


def nested_parallel(kernel_x: ConditionalLaw, kernel_y: ConditionalLaw, cost_matrix):
    assert kernel_x.T == kernel_y.T
    T = kernel_x.T
    V = [
        np.zeros([kernel_x.nc[t], kernel_y.nc[t]]) for t in range(T)
    ]  # V_t(x_{1:t},y_{1:t})
    for t in range(T - 1, -1, -1):
        n_processes = 8 if t > 0 else 1
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

    nested_ot_value = V[0][0, 0]
    return nested_ot_value


def nested_ot_solver_py(X, Y, delta_n, markovian, parallel):
    adaptedX = path2adaptedpath(X, delta_n)
    adaptedY = path2adaptedpath(Y, delta_n)

    # Quantization map
    q2v = np.unique(np.concatenate([adaptedX, adaptedY], axis=0))
    v2q = {k: v for v, k in enumerate(q2v)}  # Value to Quantization

    # Quantized paths
    qX = np.array([[v2q[x] for x in y] for y in adaptedX])
    qY = np.array([[v2q[x] for x in y] for y in adaptedY])

    # Sort paths and transpose to (n_sample, T+1)
    qX = sort_qpath(qX.T)
    qY = sort_qpath(qY.T)

    cost_matrix = np.square(q2v[:, None] - q2v[None, :])

    kernel_x = ConditionalLaw(qX, markovian)
    kernel_y = ConditionalLaw(qY, markovian)

    start_time = time.perf_counter()
    if parallel:
        nested_ot_value = nested_parallel(kernel_x, kernel_y, cost_matrix)
    else:
        nested_ot_value = nested(kernel_x, kernel_y, cost_matrix)

    end_time = time.perf_counter()
    print("Elapsed time (Nested OT): {:.4f} seconds".format(end_time - start_time))
    print("Numerical Nested OT value: ", nested_ot_value)

    return nested_ot_value
