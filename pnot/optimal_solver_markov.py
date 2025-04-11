# MARKOVIAN VERSIONS OF THE UTILITY FUNCTIONS AND PARALLEL SOLVER

import concurrent.futures

import concurrent
import numpy as np

from pnot.solver import qpath2mu_x
import ot


def qpath2mu_x_markovian(qpath):
    """
    Quantized Path to Conditional Measure for the markovian case.

    For markovian:
       mu_x[t][x_t] = { x_{t+1}: count, ... }

    This is simply a wrapper for qpath2mu_x with markovian=True.
    """
    return qpath2mu_x(qpath, markovian=True)


def list_repr_mu_x_markovian(mu_x, q2v, quantized_value=False):
    """
    Updated representation of mu_x for the markovian case.
    In addition to the conditioning states (mu_x_c) and weights, this function
    returns:
      - v2q_x: a mapping from each state (key) to its index in the dynamic programming table,
      - mu_x_idx: for t < T, for each conditional measure (a dict mapping next state -> count),
                   an array of next state indices (looked up from v2q_x at time t+1).
    """
    # List of conditioning states at time t (the keys in mu_x[t])
    mu_x_c = [list(mu_x_t.keys()) for mu_x_t in mu_x]
    # Count of conditioning states at each t
    mu_x_cn = [len(mu_x_t) for mu_x_t in mu_x]
    # Build a mapping for each time step: state -> index (for use in V[t+1])
    v2q_x = [{k[-1]: i for i, k in enumerate(mu_x_t.keys())} for mu_x_t in mu_x]

    # For each conditional measure at time t (a dictionary over next states),
    # record the quantized values of the next states.
    if quantized_value:
        mu_x_v = [
            [np.array(list(d.keys())) for d in mu_x_t.values()] for mu_x_t in mu_x
        ]
    else:
        mu_x_v = [
            [np.array([q2v[k] for k in d.keys()]) for d in mu_x_t.values()]
            for mu_x_t in mu_x
        ]

    # Build the weights for each conditional measure
    mu_x_w0 = [[list(d.values()) for d in mu_x_t.values()] for mu_x_t in mu_x]
    mu_x_w = [[np.array(l) / sum(l) for l in mu_x_w0_t] for mu_x_w0_t in mu_x_w0]

    # For t < T-1, compute the “next‐state index” arrays using v2q_x at time t+1.
    T = len(mu_x)
    mu_x_idx = []
    for t in range(T):
        if t < T - 1:
            # For each conditional measure d at time t, convert its keys into indices via v2q_x[t+1]
            idx_list = []
            for d in mu_x[t].values():
                keys = list(d.keys())
                idx = np.array([v2q_x[t + 1][k] for k in keys])
                idx_list.append(idx)
            mu_x_idx.append(idx_list)
        else:
            mu_x_idx.append(None)  # final time step: no next indices

    # (For markovian, one could assume one next index per state, but if there are multiple
    # transitions then the above handles the general case.)

    # (The following cumulative indices are no longer used in the cost lookup.)
    mu_x_n = [[1 for _ in mu_x_t] for mu_x_t in mu_x_c]
    mu_x_cumn = [np.arange(0, count + 1) for count in mu_x_cn]

    return mu_x_c, mu_x_cn, mu_x_v, mu_x_w, mu_x_cumn, v2q_x, mu_x_idx


def nested2_parallel_markovian(
    mu_x_cn,
    mu_x_v,
    mu_x_w,
    mu_x_idx,
    nu_y_cn,
    nu_y_v,
    nu_y_w,
    nu_y_idx,
    n_processes=6,
    power=2,
):
    """
    Updated parallel nested backward induction solver for the markovian adapted OT problem.
    Instead of using cumulative slicing (mu_x_cumn), we pass the next-state index arrays (mu_x_idx, nu_y_idx)
    so that for each conditional measure the next-state cost is added correctly.
    """
    T = len(mu_x_cn)
    V = [np.zeros([mu_x_cn[t], nu_y_cn[t]]) for t in range(T)]

    for t in range(T - 1, -1, -1):
        n_proc = n_processes if t > 0 else 1
        chunks = np.array_split(range(mu_x_cn[t]), n_proc)
        args = []
        for chunk in chunks:
            # Build x_arg: for each index in the chunk, pass its next-state array from mu_x_idx[t] (if available)
            x_indices = list(chunk)
            x_v = [mu_x_v[t][i] for i in chunk]
            x_w = [mu_x_w[t][i] for i in chunk]
            if t < T - 1:
                x_next_idx = [mu_x_idx[t][i] for i in chunk]
            else:
                x_next_idx = [None for _ in chunk]
            x_arg = [x_indices, x_v, x_w, x_next_idx]

            # Build y_arg analogously
            y_indices = list(range(nu_y_cn[t]))
            y_v = [nu_y_v[t][i] for i in range(nu_y_cn[t])]
            y_w = [nu_y_w[t][i] for i in range(nu_y_cn[t])]
            if t < T - 1:
                y_next_idx = [nu_y_idx[t][i] for i in range(nu_y_cn[t])]
            else:
                y_next_idx = [None for _ in range(nu_y_cn[t])]
            y_arg = [y_indices, y_v, y_w, y_next_idx]

            Vtplus = V[t + 1] if t < T - 1 else None
            arg = (x_arg, y_arg, Vtplus, power)
            args.append(arg)

        with concurrent.futures.ProcessPoolExecutor(max_workers=n_proc) as executor:
            Vts = executor.map(chunk_process_markov, args)
        for chunk, Vt_chunk in zip(chunks, Vts):
            V[t][chunk] = Vt_chunk

    AW_2square = V[0][0, 0]
    return AW_2square


def chunk_process_markov(arg):
    """
    Updated chunk process for the markovian case.
    The x_arg and y_arg now include an extra component: the next-state index arrays.
    """
    x_arg, y_arg, Vtplus, power = arg
    # x_arg: [indices, list of mu_x_v arrays, list of mu_x_w arrays, list of x_next_idx arrays]
    # y_arg: [indices, list of nu_y_v arrays, list of nu_y_w arrays, list of y_next_idx arrays]
    indices_x, mu_x_v_list, mu_x_w_list, mu_x_next_idx_list = x_arg
    indices_y, nu_y_v_list, nu_y_w_list, nu_y_next_idx_list = y_arg
    Vt = np.zeros((len(indices_x), len(indices_y)))
    for i, (vx, wx, x_next_idx) in enumerate(
        zip(mu_x_v_list, mu_x_w_list, mu_x_next_idx_list)
    ):
        for j, (vy, wy, y_next_idx) in enumerate(
            zip(nu_y_v_list, nu_y_w_list, nu_y_next_idx_list)
        ):
            Vt[i, j] = solve_ot_markov(
                vx, wx, x_next_idx, vy, wy, y_next_idx, Vtplus, power
            )
    return Vt


def solve_ot_markov(vx, wx, x_next_idx, vy, wy, y_next_idx, Vtplus, power):
    """
    Updated OT solver.
    If Vtplus is provided (i.e. t < T-1), the cost is augmented with the corresponding future cost.
    The next-state indices (x_next_idx, y_next_idx) are used to slice Vtplus.
    """
    cost = (vx[:, None] - vy[None, :]) ** power
    if Vtplus is not None and x_next_idx is not None and y_next_idx is not None:
        # Use np.ix_ to select the submatrix corresponding to the next-state indices.
        cost += Vtplus[np.ix_(x_next_idx, y_next_idx)]
    if len(vx) == 1 or len(vy) == 1:
        res = np.dot(wx, np.dot(cost, wy))
    else:
        res = ot.lp.emd2(wx, wy, cost)
    return res
