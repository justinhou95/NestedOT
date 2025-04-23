import numpy as np
from aot_numerics.normal_ot import (
    solve_sinkhorn,
    gurobi_2d,
    solve_pot,
    solve_unbalanced,
)

from time import time
from aot_numerics.measure import (
    get_meas_for_sinkhorn,
    get_joint_prob,
    get_full_index_markov,
    get_start_next_indices,
)

"""
This file contains main functions to solve temporal causal and bicausal optimal transport problems
The three main methods that are implemented are 1) backward induction (dynamic programming principle), 
2) linear programming (LP), and 3) causal and bicausal version of Sinkhorn's algorithm

Throughout, probability measures will be represented according to a disintegration with respect to a certain graph.
(usually the graph gives either Markovian 1->2->3 or full temporal (1->2, 1->3, 2->3) structure)
a measure representation of a measure mu on R^T corresponding to a graph g on T nodes is defined as follows:
a normal d-dimensional (marginal) probability measure is defined as a tuple of:
([n, d] np-array with x-values of the points, and a list like object with n entries p_i>= 0 which sum to one).
(mostly d=1 is used for marginals)
For top nodes (i.e., all nodes without parents), only their marginal distribution is given. These
are encoded as mu(node, []), where node is the specifier of the graph for this node.
More generally, mu(node, x_parents) returns a marginal probability measure, which is the conditional
distribution of mu given that its parents take the values x_parents.
"""


# FUNCTIONS TO SOLVE BICAUSAL PROBLEMS WITH BACKWARD INDUCTION


def make_new_vf(p_mu, p_nu, vals):
    # program used to create a new value function out of the solutions to the DPP problems at a given node
    n_vf = len(p_nu)

    def new_vf(x, y):
        ind1 = np.where((np.abs(p_mu - x) <= 10**-6).all(axis=1))
        ind1 = ind1[0][0]
        ind2 = np.where((np.abs(p_nu - y) <= 10**-6).all(axis=1))
        ind2 = ind2[0][0]
        return vals[ind1 * n_vf + ind2]

    return new_vf


def solve_dynamic(
    cost,
    mu,
    nu,
    supp_mu,
    supp_nu,
    g,
    index_mu=0,
    index_nu=0,
    outputflag=1,
    method="gurobi",
):
    """
    should solve a (graph) causal Wasserstein model using dynamic programming
    :param cost: a list of functions. each function contains two elements. first element is a list of nodes that this
    entry depends on. Second is a function that takes as input a tuple of support points of the relevant nodes and
    returns a value. The "true" cost function of the OT problem is the sum of all individual costs.
    :param mu: first marginal. see DAGmeasures.py regarding structure
    :param nu: same as mu
    :param supp_mu: support of first marginal. See DAGmeasures.py regarding structure
    :param supp_nu: same as supp_mu
    :param g: the graph structure
    :param method: 'gurobi', 'pot', 'sinkhorn', 'unbalanced' specifies which method is used to solve OT problems
    :return:
    """
    out_vals = []
    ordering = g.topologicalSort()
    T = g.V
    cur_V = cost  # current list of relevant value functions
    # each element of cur_V is a list with two entries. The first entry is a list of nodes that this value function
    # depends on, and the second entry is a function which returns a value for each tuple of support points of the
    # respective nodes

    optis_nodes = []  # list with T entries
    rel_node_list = (
        []
    )  # list with T entries, each being a tuple of the relevant nodes for the current node
    optis_pars = []  # list with T entries, each being an [n_1, 2*d] array of points
    optis_meas = []  # list with T entries, each being a list of n_1 entries,
    # each being a list with two entries, one [n_2, 2] array of points and one list of length
    # n_2 with weights

    for i_t in range(T):
        # get current node that we work backwards from
        v_here = ordering[T - 1 - i_t]
        optis_nodes.append(v_here)

        # get relevant nodes that are in connection with current node
        rel_nodes = g.parents[v_here].copy()

        # get relevant nodes and value functions that are in connection with current node
        rel_value_funs = []
        cur_V_old = cur_V.copy()
        cur_V = []
        for v_old in cur_V_old:
            if v_here in v_old[0]:
                rel_nodes.extend(v_old[0])
                rel_value_funs.append(v_old)
            else:
                cur_V.append(v_old)

        # relevant nodes should only include all connected nodes, but not the current node, sorted
        rel_nodes = list(set(rel_nodes))
        if v_here in rel_nodes:
            rel_nodes.remove(v_here)
        n_rel = len(rel_nodes)
        rel_nodes.sort()
        rel_node_list.append(tuple(rel_nodes))

        # # get all parent indices ... I don't think I need this after all.
        # par_indices = []
        # for v_par in g.parents[v_here]:
        #     par_indices.append(rel_nodes.index(v_par))
        # par_indices = np.array(par_indices)
        par_indices = g.parents[v_here]

        # make another array of all relevant nodes including the current node, sorted
        rel_plus_cur_nodes = rel_nodes.copy()
        rel_plus_cur_nodes.append(v_here)
        rel_plus_cur_nodes.sort()

        # get relevant supports of all measures involved:
        supp_h_mu = supp_mu(rel_nodes)  # array of size N_h_mu x n_rel
        supp_h_nu = supp_nu(rel_nodes)  # array of size N_h_nu x n_rel
        N_h_mu = len(supp_h_mu)
        N_h_nu = len(supp_h_nu)

        # if the second dimension is for some reason dimension zero, this basically still means that there are no relevant parents
        if np.prod(np.array(supp_h_mu).shape) == 0:
            N_h_mu = 0
        if np.prod(np.array(supp_h_nu).shape) == 0:
            N_h_nu = 0

        # iterate over pairs of support points of the related nodes
        vals = []
        optis = []

        # look at the case where both mu and nu have no relevant parents
        if N_h_nu == N_h_mu == 0:
            # for each pair of support points, get the disintegrations (i.e. measures) that are relevant for the
            # OT problem
            if index_mu == 1:
                input_mu = mu(v_here, 0)
            else:
                input_mu = mu(v_here, [])
            if index_nu == 1:
                input_nu = nu(v_here, 0)
            else:
                input_nu = nu(v_here, [])

            # for each pair of support points, build the cost function out of the rel_value_funs for the OT problem
            def input_fun(x, y):
                out = 0
                for vf in rel_value_funs:
                    out += vf[1](x[0:1], y[0:1])
                return out

            # solve OT problem!
            ov, opti = solve_pot(
                            [input_mu, input_nu],
                            input_fun,
                            r_opti=1,
                            outputflag=outputflag,
                        )
            out_vals.append(ov)
            vals.append(ov)

            # set optimizer:
            optis_par = np.zeros([1, 0])  # no parent values
            pmu = input_mu[0]
            pnu = input_nu[0]
            nmu = len(pmu)
            nnu = len(pnu)
            optis_x = np.zeros([nmu * nnu, 2])
            optis_w = np.zeros(nmu * nnu)
            for i in range(len(pmu)):
                for j in range(len(pnu)):
                    optis_x[i + nmu * j, 0] = pmu[i]
                    optis_x[i + nmu * j, 1] = pnu[j]
                    optis_w[i + nmu * j] = opti[i][j]

            optis_pars.append(optis_par)
            optis_meas.append([[optis_x, optis_w]])

        else:
            optis_par = np.zeros([N_h_nu * N_h_mu, 2 * n_rel])
            optis_meas_h = []
            for i in range(N_h_mu):
                for j in range(N_h_nu):

                    # get the x and y point that is relevant for this iteration
                    p_mu_h = supp_h_mu[i, :]
                    p_nu_h = supp_h_nu[j, :]
                    p_mu_h_ext = np.zeros(T)
                    p_nu_h_ext = np.zeros(T)
                    p_mu_h_ext[rel_nodes] = p_mu_h
                    p_nu_h_ext[rel_nodes] = p_nu_h

                    # extract parents for disintegration
                    p_par_mu_h = p_mu_h_ext[par_indices]
                    p_par_nu_h = p_nu_h_ext[par_indices]

                    # for each pair of support points, get the disintegrations (i.e. measures) that are relevant for the
                    # OT problem
                    input_mu = mu(v_here, p_par_mu_h)
                    input_nu = nu(v_here, p_par_nu_h)

                    # for each pair of support points, build cost function out of the rel_value_funs for the OT problem
                    def input_fun(x, y):
                        out = 0
                        p_mu_h_ext_vf = p_mu_h_ext.copy()
                        p_nu_h_ext_vf = p_nu_h_ext.copy()
                        p_mu_h_ext_vf[v_here] = x
                        p_nu_h_ext_vf[v_here] = y

                        for vf in rel_value_funs:
                            inds_vf = vf[0]
                            xinpvf = p_mu_h_ext_vf[inds_vf]
                            yinpvf = p_nu_h_ext_vf[inds_vf]
                            out += vf[1](xinpvf, yinpvf)
                        return out

                    # solve OT problem!
                    if method == "gurobi":
                        ov, opti = gurobi_2d(
                            [input_mu, input_nu],
                            input_fun,
                            r_opti=1,
                            outputflag=outputflag,
                        )
                    elif method == "unbalanced":
                        ov, opti = solve_unbalanced(
                            [input_mu, input_nu],
                            input_fun,
                            r_opti=1,
                            outputflag=outputflag,
                        )
                    elif method == "pot":
                        ov, opti = solve_pot(
                            [input_mu, input_nu],
                            input_fun,
                            r_opti=1,
                            outputflag=outputflag,
                        )
                    elif method == "sinkhorn":
                        ov, opti = solve_sinkhorn(
                            [input_mu, input_nu],
                            input_fun,
                            r_opti=1,
                            outputflag=outputflag,
                        )

                    vals.append(ov)
                    optis.append(opti)
                    opti = np.array(opti)

                    # set optis_x:
                    optis_par[i * N_h_nu + j, :n_rel] = p_mu_h
                    optis_par[i * N_h_nu + j, n_rel:] = p_nu_h
                    pmu = input_mu[0]
                    pnu = input_nu[0]
                    nmu = len(pmu)
                    nnu = len(pnu)
                    optis_x = np.zeros([nmu * nnu, 2])
                    optis_w = np.zeros(nmu * nnu)
                    for ii in range(len(pmu)):
                        for jj in range(len(pnu)):
                            optis_x[ii + nmu * jj, 0] = pmu[ii]
                            optis_x[ii + nmu * jj, 1] = pnu[jj]
                            optis_w[ii + nmu * jj] = opti[ii, jj]
                    optis_meas_h.append([optis_x, optis_w])
            optis_pars.append(optis_par)
            optis_meas.append(optis_meas_h)

        # build new value function out of solutions to OT problems and add it to cur_V
        new_vf = make_new_vf(supp_h_mu, supp_h_nu, vals)
        V_new = [rel_nodes, new_vf]
        cur_V.append(V_new)

        # save optimal coupling in some suitable fashion...

    # at the end, cur_V should only contain one function that does not depend on any input. The output value is the
    # returns value

    return out_vals, [optis_nodes, rel_node_list, optis_pars, optis_meas]


# FUNCTION TO DIRECTLY SOLVE CAUSAL AND BICAUSAL OT VIA LINEAR PROGRAMMING
def gurobi_bm(
    margs,
    f,
    p_dist=2,
    radial_cost=0,
    f_id=0,
    minmax="min",
    r_opti=0,
    outputflag=1,
    causal=0,
    anticausal=0,
):
    """
    :param margs: list with 2 entries, each entry being a discrete probability measure on R^n, where x_list is an [N, n] array
    :param f: function that takes two inputs, x, y, where the inputs are of the form as in the representation of the
    points in margs. Returns a single value
    :param p_dist: if radial cost is used, then this describes the Lp norm which is used.
    :param radial_cost: If 1, then f takes an arbitrary number of inputs but treats them element-wise. Each element
    which will be \|x-y\|_{p_dist} for some x, y. This allows for a faster computation of the cost matrix.
    :param f_id: if non-zero and raidal_cost nonzero, then f will be treated as the identity function.
    :param minmax: if 'min', then we minimize objective, else, we maximize
    :param r_opti: if 0, does not return optimizer. if 1, it does
    :return: optimal value (and optimizer) of the OT problem
    """
    # get relevant data from input:
    m1 = margs[0]
    m2 = margs[1]
    xl_1 = np.array(m1[0])
    xl_2 = np.array(m2[0])
    pl_1 = m1[1]
    pl_2 = m2[1]
    n1, n_dim = xl_1.shape
    n2 = len(xl_2)

    if len(xl_1.shape) == 1:
        xl_1 = xl_1.reshape(-1, 1)
    if len(xl_2.shape) == 1:
        xl_2 = xl_2.reshape(-1, 1)

    # build cost matrix:
    # print('Building cost matrix...')
    if radial_cost == 0:
        cost_mat = np.zeros([n1, n2])
        for i in range(n1):
            for j in range(n2):
                cost_mat[i, j] = f(xl_1[i, :], xl_2[j, :])
    else:
        cost_mat = np.linalg.norm(
            xl_1[:, None, :] - xl_2[None, :, :], axis=-1, ord=p_dist
        )
        if f_id == 0:
            cost_mat = f(cost_mat)

    # initialize model
    # print('Initializing model...')
    m = Model("Primal")
    if outputflag == 0:
        m.setParam("OutputFlag", 0)
    pi_var = m.addVars(n1, n2, lb=0, ub=1, name="pi_var")

    # add marginal constraints
    # print('Adding constraints...')
    m.addConstrs((pi_var.sum(i, "*") == pl_1[i] for i in range(n1)), name="first_marg")
    m.addConstrs((pi_var.sum("*", i) == pl_2[i] for i in range(n2)), name="second_marg")

    # add causal constraint: (Note: doesn't seem very efficient, but not sure how else to do)
    causal_count = 0
    if causal == 1:
        for t in range(1, n_dim):
            x_t_arr, ind_inv = np.unique(xl_1[:, :t], axis=0, return_inverse=True)
            for ind_t in range(len(x_t_arr)):
                pos_h = np.where(ind_inv == ind_t)[0]
                y_t_arr, ind_inv_y = np.unique(xl_2[:, :t], axis=0, return_inverse=True)
                for ind_t_y in range(len(y_t_arr)):
                    pos_h_y = np.where(ind_inv_y == ind_t_y)[0]
                    x_tp_arr, ind_inv_p = np.unique(
                        xl_1[pos_h, : t + 1], axis=0, return_inverse=True
                    )
                    for ind_xp in range(len(x_tp_arr)):
                        pos_xtp = np.where(ind_inv_p == ind_xp)[0]
                        pos_xtp_real = pos_h[pos_xtp]
                        pi_sum_left = 0
                        for i_x in pos_xtp_real:
                            for i_y in pos_h_y:
                                pi_sum_left += pi_var[i_x, i_y]
                        pi_sum_right = 0
                        for i_x in pos_h:
                            for i_y in pos_h_y:
                                pi_sum_right += pi_var[i_x, i_y]
                        mu_sum_left = 0
                        for i_x in pos_h:
                            mu_sum_left += pl_1[i_x]
                        mu_sum_right = 0
                        for i_x in pos_xtp_real:
                            mu_sum_right += pl_1[i_x]

                        causal_count += 1
                        m.addConstr(
                            pi_sum_left * mu_sum_left == pi_sum_right * mu_sum_right,
                            name="causal_"
                            + str(t)
                            + "_"
                            + str(ind_t)
                            + "_"
                            + str(ind_t_y)
                            + "_"
                            + str(ind_xp),
                        )

    if anticausal == 1:
        for t in range(1, n_dim):
            x_t_arr, ind_inv = np.unique(xl_2[:, :t], axis=0, return_inverse=True)
            for ind_t in range(len(x_t_arr)):
                pos_h = np.where(ind_inv == ind_t)[0]

                y_t_arr, ind_inv_y = np.unique(xl_1[:, :t], axis=0, return_inverse=True)
                for ind_t_y in range(len(y_t_arr)):
                    pos_h_y = np.where(ind_inv_y == ind_t_y)[0]

                    x_tp_arr, ind_inv_p = np.unique(
                        xl_2[pos_h, : t + 1], axis=0, return_inverse=True
                    )
                    # TODO: note that we have to concatenate pos_h and pos_p to get real index! (done, but good to keep in mind)

                    for ind_xp in range(len(x_tp_arr)):
                        pos_xtp = np.where(ind_inv_p == ind_xp)[0]
                        pos_xtp_real = pos_h[pos_xtp]

                        pi_sum_left = 0
                        for i_x in pos_xtp_real:
                            for i_y in pos_h_y:
                                pi_sum_left += pi_var[i_y, i_x]

                        pi_sum_right = 0
                        for i_x in pos_h:
                            for i_y in pos_h_y:
                                pi_sum_right += pi_var[i_y, i_x]

                        mu_sum_left = 0
                        for i_x in pos_h:
                            mu_sum_left += pl_2[i_x]

                        mu_sum_right = 0
                        for i_x in pos_xtp_real:
                            mu_sum_right += pl_2[i_x]

                        m.addConstr(
                            pi_sum_left * mu_sum_left == pi_sum_right * mu_sum_right,
                            name="anticausal_"
                            + str(t)
                            + "_"
                            + str(ind_t)
                            + "_"
                            + str(ind_t_y)
                            + "_"
                            + str(ind_xp),
                        )

    # Specify objective function
    if minmax == "min":
        obj = quicksum(
            [cost_mat[i, j] * pi_var[i, j] for i in range(n1) for j in range(n2)]
        )
        m.setObjective(obj, GRB.MINIMIZE)
    else:
        obj = quicksum(
            [cost_mat[i, j] * pi_var[i, j] for i in range(n1) for j in range(n2)]
        )
        m.setObjective(obj, GRB.MAXIMIZE)

    # solve model
    m.optimize()
    objective_val = m.ObjVal

    if r_opti == 0:
        return objective_val
    else:
        return objective_val, [[pi_var[i, j].x for j in range(n2)] for i in range(n1)]


def sinkhorn_bicausal_markov(
    mu_list,
    nu_list,
    cost_list,
    n_list,
    m_list,
    eps_stop=10**-4,
    max_iter=10**4,
    reshape=True,
    outputflag=0,
):
    # Only for MARKOV - MARKOV marginals, bicausal!
    """

    :param mu_list: as output by get_meas_for_sinkhorn
    :param nu_list: as output by get_meas_for_sinkhorn
    :param cost_list: list of matrices, one for each time point (markov case). Notably, the cost functions should
                    already be kernelized, i.e., values are exp(-c) instead of c
    :param n_list: sizes of supports for mu for each time step
    :param m_list: sizes of supports for nu for each time step
    :return:
    """
    t_max = len(mu_list)

    # initializing dual functions. We specify them in a multiplicative way, i.e., compared to the paper, we store values
    # of exp(f_t) and exp(g_t) instead of f_t and g_t, which is in line with standard implementations of Sinkhorn's
    tinit = time()
    f_1 = np.ones(n_list[0])
    g_1 = np.ones(m_list[0])
    f_list = [f_1]
    g_list = [g_1]
    const_f_list = [0]
    const_g_list = [0]
    for t in range(1, t_max):
        f_h = [
            [np.ones([len(mu_list[t][i][1]), 1]) for j in range(m_list[t - 1])]
            for i in range(n_list[t - 1])
        ]
        g_h = [
            [np.ones([1, len(nu_list[t][j][1])]) for j in range(m_list[t - 1])]
            for i in range(n_list[t - 1])
        ]
        c_f_h = [[1 for j in range(m_list[t - 1])] for i in range(n_list[t - 1])]
        c_g_h = [[1 for j in range(m_list[t - 1])] for i in range(n_list[t - 1])]
        f_list.append(f_h)
        g_list.append(g_h)
        const_f_list.append(c_f_h)
        const_g_list.append(c_g_h)
    if outputflag:
        print("Initializing took " + str(time() - tinit) + " seconds")

    # Define update iterations:
    t_funs = time()

    def update_f_t(mut, nut, gt, ct):
        """

        :param mut: should be of shape (a, 1)
        :param nut: should be of shape (1, b)
        :param gt: should be of shape (1, b)
        :param ct: should be of shape (a, b)
        :return: array of shape (a, 1) representing f_t
        """
        # at = 1. / np.sum(gt * ct * nut, axis=1, keepdims=True)
        # at = 1. / np.dot(ct, (gt*nut).T)
        at = 1.0 / np.matmul(ct, (gt * nut).T)
        cth = np.sum(np.log(at) * mut)
        return at / np.exp(cth), cth

    def update_g_t(mut, nut, ft, ct):
        """

        :param mut: should be of shape (a, 1)
        :param nut: should be of shape (1, b)
        :param ft: should be of shape (a, 1)
        :param ct: should be of shape (a, b)
        :return: array of shape (1, b) representing g_t
        """
        # bt = 1. / np.sum(ft*ct*mut, axis=0, keepdims=True)
        # bt = 1. / np.dot(ct.T, ft*mut).T
        bt = 1.0 / np.matmul(ct.T, ft * mut).T
        cth = np.sum(np.log(bt) * nut)
        return bt / np.exp(cth), cth

    def update_f_1(mut, nut, gt, ct):
        # inputs as for update_f_t
        at = 1.0 / np.sum(gt * ct * nut, axis=1, keepdims=True)
        return at, np.sum(np.log(at) * mut)

    def update_g_1(mut, nut, ft, ct):
        # inputs as for update_g_t
        bt = 1.0 / np.sum(ft * ct * mut, axis=0, keepdims=True)
        return bt, np.sum(np.log(bt) * nut)

    def full_update_f_list():
        for t_m in range(t_max):
            t = t_max - t_m - 1
            if t > 0:
                cvnew = np.ones([n_list[t - 1], m_list[t - 1]])
            if t == 0:
                f_list[0], value_f = update_f_1(
                    mu_list[0][1],
                    nu_list[0][1],
                    g_list[0],
                    cost_list[0][mu_list[0][0], :][:, nu_list[0][0]]
                    * cvh[mu_list[0][0], :][:, nu_list[0][0]],
                )
            elif t == t_max - 1:
                for i in range(n_list[t - 1]):
                    for j in range(m_list[t - 1]):
                        f_list[t][i][j], cvnew[i, j] = update_f_t(
                            mu_list[t][i][1],
                            nu_list[t][j][1],
                            g_list[t][i][j],
                            cost_list[t][mu_list[t][i][0], :][:, nu_list[t][j][0]],
                        )
            else:
                for i in range(n_list[t - 1]):
                    for j in range(m_list[t - 1]):
                        f_list[t][i][j], cvnew[i, j] = update_f_t(
                            mu_list[t][i][1],
                            nu_list[t][j][1],
                            g_list[t][i][j],
                            cost_list[t][mu_list[t][i][0], :][:, nu_list[t][j][0]]
                            * cvh[mu_list[t][i][0], :][:, nu_list[t][j][0]],
                        )
            cvh = np.exp(-cvnew.copy())
            const_f_list[t] = cvh.copy()
        return value_f

    def full_update_g_list():
        for t_m in range(t_max):
            t = t_max - t_m - 1
            if t > 0:
                cvnew = np.ones([n_list[t - 1], m_list[t - 1]])
            if t == 0:
                g_list[0], value_g = update_g_1(
                    mu_list[0][1],
                    nu_list[0][1],
                    f_list[0],
                    cost_list[0][mu_list[0][0], :][:, nu_list[0][0]]
                    * cvh[mu_list[0][0], :][:, nu_list[0][0]],
                )
            elif t == t_max - 1:
                for i in range(n_list[t - 1]):
                    for j in range(m_list[t - 1]):
                        g_list[t][i][j], cvnew[i, j] = update_g_t(
                            mu_list[t][i][1],
                            nu_list[t][j][1],
                            f_list[t][i][j],
                            cost_list[t][mu_list[t][i][0], :][:, nu_list[t][j][0]],
                        )
            else:
                for i in range(n_list[t - 1]):
                    for j in range(m_list[t - 1]):
                        g_list[t][i][j], cvnew[i, j] = update_g_t(
                            mu_list[t][i][1],
                            nu_list[t][j][1],
                            f_list[t][i][j],
                            cost_list[t][mu_list[t][i][0], :][:, nu_list[t][j][0]]
                            * cvh[mu_list[t][i][0], :][:, nu_list[t][j][0]],
                        )
            cvh = np.exp(-cvnew.copy())
            const_g_list[t] = cvh.copy()
        return value_g

    if outputflag:
        print("Defining update functions took " + str(time() - t_funs) + " seconds")

    if reshape:
        # reshape inputs
        # we want mu_list[t][i][1] to be shaped (a, 1) and nu_list[t][j][1] to be shaped (1, b) for some a and b that may
        # depend on i and j
        t_reshape = time()
        for t in range(t_max):
            if t == 0:
                if len(mu_list[t][1].shape) == 1:
                    mu_list[t][1] = np.expand_dims(mu_list[t][1], 1)
                if len(nu_list[t][1].shape) == 1:
                    nu_list[t][1] = np.expand_dims(nu_list[t][1], 0)
                if len(mu_list[t]) == 2:
                    mu_list[t].append(np.log(mu_list[t][1]))
                if len(nu_list[t]) == 2:
                    nu_list[t].append(np.log(nu_list[t][1]))
            else:
                for i in range(n_list[t - 1]):
                    if len(mu_list[t][i][1].shape) == 1:
                        mu_list[t][i][1] = np.expand_dims(mu_list[t][i][1], 1)
                    if len(mu_list[t][i]) == 2:
                        mu_list[t][i].append(np.log(mu_list[t][i][1]))

                for j in range(m_list[t - 1]):
                    if len(nu_list[t][j][1].shape) == 1:
                        nu_list[t][j][1] = np.expand_dims(nu_list[t][j][1], 0)
                    if len(nu_list[t][j]) == 2:
                        nu_list[t][j].append(np.log(nu_list[t][j][1]))

        if outputflag:
            print("Reshaping input took " + str(time() - t_reshape) + " seconds")

    t_solve = time()
    prev_val = -(10**8)
    value_f = -100
    value_g = -100
    iter_h = 0
    while iter_h < max_iter and np.abs(prev_val - value_f - value_g) > eps_stop:
        if iter_h % 10 == 0 and outputflag:
            print(
                "Current iteration:",
                iter_h,
                "Current value:",
                value_f + value_g,
                "Current time:",
                time() - t_solve,
            )
        iter_h += 1
        prev_val = value_f + value_g
        value_f = full_update_f_list()
        value_g = full_update_g_list()
        # print(value_f)
        # print(value_g)
    if outputflag:
        print("Solving took " + str(time() - t_solve) + " seconds")

    # get value without entropy
    for t_m in range(t_max):
        t = t_max - t_m - 1
        if t > 0:
            V_t = np.zeros([n_list[t - 1], m_list[t - 1]])
        if t == t_max - 1:
            for i in range(n_list[t - 1]):
                for j in range(m_list[t - 1]):
                    V_t[i, j] = np.sum(
                        -np.log(cost_list[t][mu_list[t][i][0], :][:, nu_list[t][j][0]])
                        * f_list[t][i][j]
                        * g_list[t][i][j]
                        * cost_list[t][mu_list[t][i][0], :][:, nu_list[t][j][0]]
                        * (1.0 / const_g_list[t][i, j])
                        * mu_list[t][i][1]
                        * nu_list[t][j][1]
                    )
        elif t > 0:
            for i in range(n_list[t - 1]):
                for j in range(m_list[t - 1]):
                    V_t[i, j] = np.sum(
                        (
                            -np.log(
                                cost_list[t][mu_list[t][i][0], :][:, nu_list[t][j][0]]
                            )
                            + V_tp[mu_list[t][i][0], :][:, nu_list[t][j][0]]
                        )
                        * f_list[t][i][j]
                        * g_list[t][i][j]
                        * cost_list[t][mu_list[t][i][0], :][:, nu_list[t][j][0]]
                        * const_g_list[t + 1][mu_list[t][i][0], :][:, nu_list[t][j][0]]
                        * (1.0 / const_g_list[t][i, j])
                        * mu_list[t][i][1]
                        * nu_list[t][j][1]
                    )
        else:
            value = np.sum(
                (
                    -np.log(cost_list[0][mu_list[0][0], :][:, nu_list[0][0]])
                    + V_tp[mu_list[t][0], :][:, nu_list[t][0]]
                )
                * f_list[0]
                * g_list[0]
                * cost_list[0][mu_list[0][0], :][:, nu_list[0][0]]
                * const_g_list[t + 1][mu_list[t][0], :][:, nu_list[t][0]]
                * mu_list[0][1]
                * nu_list[0][1]
            )
        V_tp = V_t.copy()
    return value
    # return 0.0


def sinkhorn_causal_markov(
    mu_list,
    nu_list,
    cost_list,
    n_list,
    m_list,
    nu_index_full,
    nu_next_list,
    nu_probs_tm1,
    eps_stop=10**-4,
    max_iter=10**4,
    outputflag=0,
    reshape=1,
):
    # For MARKOV - MARKOV marginals, causal!
    """

    :param mu_list: as output by get_meas_for_sinkhorn
    :param nu_list: as output by get_meas_for_sinkhorn
    :param cost_list: list of matrices, one for each time point (markov case). Already saved in kernel style, i.e.,
                        exp(-c) is given instead of c
    :param n_list: sizes of supports for mu for each time step
    :param m_list: sizes of supports for nu for each time step
    :return:
    """
    t_max = len(mu_list)

    # initializing functions
    tinit = time()
    f_1 = np.ones(n_list[0])
    f_list = [f_1]
    for t in range(1, t_max):
        f_h = [
            [np.ones([len(mu_list[t][i][1]), 1]) for j in nu_index_full[t - 1]]
            for i in range(n_list[t - 1])
        ]
        f_list.append(f_h)
    s_nu_1Tm1 = nu_index_full[-2]
    g = [np.ones([1, len(nu_list[t_max - 1][iTm1[-2]][1])]) for iTm1 in s_nu_1Tm1]
    g_vals = [0 for i in range(t_max)]
    if outputflag:
        print("Initializing took " + str(time() - tinit) + " seconds")

    # Define update iterations:
    t_funs = time()

    def update_f_capt(mut, nut, gt, ct):
        # step for t_max, where gt = g
        # input is conditional on index of x_{T-1}, y_{1:T-1}
        # gt should be shaped (1, b), ct (a, b), mut (a, 1), nut (1, b)

        # at = 1. / np.sum(gt * ct * nut, axis=1, keepdims=True)
        at = 1.0 / np.dot(ct, (gt * nut).T)
        cth = np.sum(np.log(at) * mut)
        return at / np.exp(cth), cth

    def update_f_t(mut, nut, vtp, ct):
        """

        :param mut: should be of shape (a, 1)
        :param nut: should be of shape (1, b)
        :param gt: should be of shape (1, b)
        :param ct: should be of shape (a, b)
        :return: array of shape (a, 1) representing f_t
        """
        # at = 1. / np.sum(vtp * ct * nut, axis=1, keepdims=True)
        at = 1.0 / np.dot(ct * vtp, nut.T)

        cth = np.sum(np.log(at) * mut)
        return at / np.exp(cth), cth

    def update_f_1(mut, nut, vtp, ct):
        # inputs as for update_f_t
        at = 1.0 / np.sum(vtp * ct * nut, axis=1, keepdims=True)
        return at, np.sum(np.log(at) * mut)

    def full_update_f_list():
        for t_m in range(t_max):
            t = t_max - t_m - 1
            if t > 0:
                cvnew = np.zeros([n_list[t - 1], len(nu_index_full[t - 1])])
            if t == 0:
                f_list[0], value_f = update_f_1(
                    mu_list[0][1],
                    nu_list[0][1],
                    cvh,
                    cost_list[0][mu_list[0][0], :][:, nu_list[0][0]],
                )
            elif t == t_max - 1:
                for i in range(n_list[t - 1]):
                    for j, ind_full in enumerate(nu_index_full[t - 1]):
                        f_list[t][i][j], cvnew[i, j] = update_f_capt(
                            mu_list[t][i][1],
                            nu_list[t][ind_full[-2]][1],
                            g[j],
                            cost_list[t][mu_list[t][i][0], :][
                                :, nu_list[t][ind_full[-2]][0]
                            ],
                        )
            else:
                for i in range(n_list[t - 1]):
                    for j, ind_full in enumerate(nu_index_full[t - 1]):
                        isn_start = nu_next_list[t - 1][j]
                        isn_stop = nu_next_list[t - 1][j + 1]
                        f_list[t][i][j], cvnew[i, j] = update_f_t(
                            mu_list[t][i][1],
                            nu_list[t][ind_full[-2]][1],
                            cvh[mu_list[t][i][0], :][:, isn_start:isn_stop],
                            cost_list[t][mu_list[t][i][0], :][
                                :, nu_list[t][ind_full[-2]][0]
                            ],
                        )
            cvh = np.exp(-cvnew.copy())
        return value_f

    def get_vnext(vcur, ft, ct, mut):
        """
        input is conditional on x_{t-1} and y_{1:T-1}
        :param vcur: of shape (a, b) (for t=T, simply a scalar with value 1)
        :param ft: of shape (a, b) if t=T, otherwise of shape (a, 1)
        :param ct: of shape (a, b) if t=T, otherwise of shape (a, 1)
        :param mut: of shape (a, 1)
        :return: vector of shape (b) representing one element of vcur for the next step
        """
        # return np.sum(vcur * ft * ct * mut, axis=0)
        return np.dot((vcur * ct).T, (ft * mut)).flatten()

    def update_g():
        itm1 = nu_index_full[-2]
        val_g = 0
        for (
            ind_tot
        ) in (
            itm1
        ):  # iterate over y_1, ..., y_{T-1} jointly. So g(y) = g(y_{1:T-1}, y_T)
            snu_shape = np.prod(np.shape(nu_list[t_max - 1][ind_tot[-2]][1]))
            Vcur = 1
            for t_m in range(t_max):
                t = t_max - t_m - 1
                if t > 1:
                    smu_shape = n_list[t - 1]
                else:
                    smu_shape = 1
                Vnext = np.zeros([smu_shape, snu_shape])

                for i in range(smu_shape):
                    if t == t_max - 1:
                        Vnext[i, :] = get_vnext(
                            Vcur,
                            f_list[t][i][ind_tot[-1]],
                            cost_list[t][mu_list[t][i][0], :][
                                :, nu_list[t][ind_tot[-2]][0]
                            ],
                            mu_list[t][i][1],
                        )
                    elif t > 0:
                        Vnext[i, :] = get_vnext(
                            Vcur[mu_list[t][i][0], :],
                            f_list[t][i][ind_tot[2 * t - 1]],
                            cost_list[t][
                                mu_list[t][i][0], ind_tot[2 * t] : ind_tot[2 * t] + 1
                            ],
                            mu_list[t][i][1],
                        )
                    else:  # t = 0
                        Vnext[i, :] = get_vnext(
                            Vcur[0:1, :],
                            f_list[0],
                            cost_list[t][
                                mu_list[t][i][0], ind_tot[2 * t] : ind_tot[2 * t] + 1
                            ],
                            mu_list[t][1],
                        )

                Vcur = Vnext.copy()
                g_vals[t] = Vcur.copy()
            g[ind_tot[-1]] = 1.0 / Vcur.copy()
            val_g += nu_probs_tm1[ind_tot[-1]] * np.sum(
                np.log(g[ind_tot[-1]]) * nu_list[t_max - 1][ind_tot[-2]][1]
            )
        return val_g

    if outputflag:
        print("Defining update functions took " + str(time() - t_funs) + " seconds")

    # reshape inputs
    # we want mu_list[t][i][1] to be shaped (a, 1) and nu_list[t][j][1] to be shaped (1, b) for some a and b that may
    # depend on i and j
    if reshape:
        t_reshape = time()
        for t in range(t_max):
            if t == 0:
                mu_list[t][1] = np.expand_dims(mu_list[t][1], 1)
                nu_list[t][1] = np.expand_dims(nu_list[t][1], 0)
            else:
                for i in range(n_list[t - 1]):
                    mu_list[t][i][1] = np.expand_dims(mu_list[t][i][1], 1)
                for j in range(m_list[t - 1]):
                    nu_list[t][j][1] = np.expand_dims(nu_list[t][j][1], 0)
        if outputflag:
            print("Reshaping input took " + str(time() - t_reshape) + " seconds")

    t_solve = time()
    prev_val = -(10**8)
    value_f = -100
    value_g = -100
    iter_h = 0
    while iter_h < max_iter and np.abs(prev_val - value_f - value_g) > eps_stop:
        if iter_h % 10 == 0 and outputflag:
            print("Current iteration:", iter_h, "Current value:", value_f + value_g)
        iter_h += 1
        prev_val = value_f + value_g
        value_f = full_update_f_list()
        value_g = update_g()
        # print(value_f)
        # print(value_g)
    if outputflag:
        print("Solving took " + str(time() - t_solve) + " seconds")

    # get value without entropy
    itm1 = nu_index_full[-2]
    value = 0
    const_val = 0
    for j, ind_tot in enumerate(
        itm1
    ):  # iterate over y_1, ..., y_{T-1} jointly. So g(y) = g(y_{1:T-1}, y_T)
        snu_shape = np.prod(np.shape(nu_list[t_max - 1][ind_tot[-2]][1]))
        for t_m in range(t_max):
            t = t_max - t_m - 1
            if t > 1:
                smu_shape = n_list[t - 1]
            else:
                smu_shape = 1
            Vnext = np.zeros([smu_shape, snu_shape])
            Vconst = np.zeros([smu_shape, snu_shape])
            for i in range(smu_shape):
                if t == t_max - 1:
                    Vconst[i, :] = np.sum(
                        cost_list[t][mu_list[t][i][0], :][:, nu_list[t][ind_tot[-2]][0]]
                        * f_list[t][i][ind_tot[-1]]
                        * mu_list[t][i][1],
                        axis=0,
                    )
                    Vnext[i, :] = np.sum(
                        -np.log(
                            cost_list[t][mu_list[t][i][0], :][
                                :, nu_list[t][ind_tot[-2]][0]
                            ]
                        )
                        * cost_list[t][mu_list[t][i][0], :][
                            :, nu_list[t][ind_tot[-2]][0]
                        ]
                        * f_list[t][i][ind_tot[-1]]
                        * mu_list[t][i][1],
                        axis=0,
                    )
                elif t > 0:
                    Vconst[i, :] = np.sum(
                        Vconst_cur[mu_list[t][i][0], :]
                        * cost_list[t][
                            mu_list[t][i][0], ind_tot[2 * t] : ind_tot[2 * t] + 1
                        ]
                        * f_list[t][i][ind_tot[2 * t - 1]]
                        * mu_list[t][i][1],
                        axis=0,
                    )
                    Vnext[i, :] = np.sum(
                        (
                            Vcur[mu_list[t][i][0], :]
                            + Vconst_cur[mu_list[t][i][0], :]
                            * (
                                -np.log(
                                    cost_list[t][
                                        mu_list[t][i][0],
                                        ind_tot[2 * t] : ind_tot[2 * t] + 1,
                                    ]
                                )
                            )
                        )
                        * cost_list[t][
                            mu_list[t][i][0], ind_tot[2 * t] : ind_tot[2 * t] + 1
                        ]
                        * f_list[t][i][ind_tot[2 * t - 1]]
                        * mu_list[t][i][1],
                        axis=0,
                    )
                else:  # t = 0
                    Vconst[i, :] = np.sum(
                        Vconst_cur[mu_list[t][0], :]
                        * cost_list[t][
                            mu_list[t][0], ind_tot[2 * t] : ind_tot[2 * t] + 1
                        ]
                        * f_list[0]
                        * mu_list[t][1],
                        axis=0,
                    )
                    Vnext[i, :] = np.sum(
                        (
                            Vcur[mu_list[t][0], :]
                            + Vconst_cur[mu_list[t][i][0], :]
                            * (
                                -np.log(
                                    cost_list[t][
                                        mu_list[t][0],
                                        ind_tot[2 * t] : ind_tot[2 * t] + 1,
                                    ]
                                )
                            )
                        )
                        * cost_list[t][
                            mu_list[t][0], ind_tot[2 * t] : ind_tot[2 * t] + 1
                        ]
                        * f_list[0]
                        * mu_list[t][1],
                        axis=0,
                    )
            Vcur = Vnext.copy()
            Vconst_cur = Vconst.copy()
        value += nu_probs_tm1[ind_tot[-1]] * np.sum(
            Vcur * nu_list[t_max - 1][ind_tot[-2]][1] * g[j]
        )
        const_val += nu_probs_tm1[ind_tot[-1]] * np.sum(
            Vconst_cur * nu_list[t_max - 1][ind_tot[-2]][1] * g[j]
        )

    return value
