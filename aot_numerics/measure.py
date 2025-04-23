import numpy as np
from sklearn.cluster import KMeans
from itertools import product
from collections import defaultdict


"""
This file contains functions to obtain different representations of the same probability measure, to generate graphs, 
to generate some simple probability measures, and to obtain an adapted version of an empirical measure.
"""


# FUNCTIONS TO RESHAPE MEASURES FOR CAUSAL AND BICAUSAL VERSION OF SINKHORN'S ALGORITHM
def get_meas_for_sinkhorn(mu, supp_mu, t_max, markov=1):
    """
    Only implemented for Markovian measures at the moment
    mu and mu_supp are as usual for DPP
    Returns:
    x_list gives the supports at each time step
    a simply measure will be represented by a list of indices, and an array of weights, corresponding to the indices
    so mu_list[0] is just such a simply measure
    On the other hand, mu_list[i] for i > 0 will be a list, each entry k corresponding to the mu_i(x_k, \cdot),
    so mu_list[i] is basically the kernel mu_i(\cdot, \cdot), and each entry k specifies the value (i.e., the measure)
    """

    x_list = [supp_mu([i]).flatten() for i in range(t_max)]
    mu_list = []
    mu0_s_pre = mu(0, [])[0].flatten()
    w0_pre = np.array(mu(0, [])[1])
    mu0 = []
    w0 = []
    for ind in range(len(mu0_s_pre)):
        for ind2, x_h in enumerate(x_list[0]):
            if mu0_s_pre[ind] == x_h:
                if ind2 in mu0:
                    i_where = mu0.index(ind2)
                    w0[i_where] += w0_pre[ind]
                else:
                    mu0.append(ind2)
                    w0.append(w0_pre[ind])
                break
    mu_list.append([tuple(mu0), np.array(w0)])

    for t in range(1, t_max):
        mu_h = []
        for ind_minus, x_minus in enumerate(x_list[t - 1]):
            mut_s_pre = mu(t, [x_minus])[0].flatten()
            wt_pre = np.array(mu(t, [x_minus])[1])
            mut = []
            wt = []
            for ind in range(len(mut_s_pre)):
                for ind2, x_h in enumerate(x_list[t]):
                    if mut_s_pre[ind] == x_h:
                        if ind2 in mut:
                            i_where = mut.index(ind2)
                            wt[i_where] += wt_pre[ind]
                        else:
                            mut.append(ind2)
                            wt.append(wt_pre[ind])
            mu_h.append([tuple(mut), np.array(wt)])
        mu_list.append(mu_h)

    return x_list, mu_list


def get_full_index_markov(mu_list):
    """
    should return a list of the same length as mu_list, where each entry represents the support of mu_{1:t} in the
    "full-index"-version. I.e., each entry of out[t], i.e., out[t][k], is a list where the j=0 entry is ind_0,
    the j=1 entry is ind_{0:0}, the j=2 entry is ind_1, the j=3 entry is ind_{0:1}, ...,
    the 2j entry is ind_{j} and the 2j+1 entry is ind_{0:j} (for j up to t). k ranges through the number of different
    support points of mu_{1:t}
    """
    if len(mu_list) == 1:
        return [
            [[mu_list[0][0][i], mu_list[0][0][i]] for i in range(len(mu_list[0][0]))]
        ]
    else:
        glprev_full = get_full_index_markov(mu_list[:-1])
        glprev = glprev_full[-1]
        out_l = []
        itot_1t = 0
        for i_tot_prev in glprev:
            i_sup_next = i_tot_prev[-2]
            sh_list = mu_list[-1][i_sup_next][0]
            for sh in sh_list:
                out_l_add = i_tot_prev.copy()
                out_l_add.append(sh)
                out_l_add.append(itot_1t)
                itot_1t += 1
                out_l.append(out_l_add.copy())
        gl_out = glprev_full.copy()
        gl_out.append(out_l)
        return gl_out


def get_start_next_indices(gl_out):
    # for causal version of Sinkhorn's algorithm.
    # output is a list with t_max-1 entries, each being a list of numbers, say out_t
    # out_t[i]:out_t[i+1] gives the range of indices for the joint support of mu_{1:t+1} which share the same elements
    # x_{1:t} (so in this range of indices, only x_{t+1} varies).
    t_max = len(gl_out)
    out_next_l = []
    for t in range(1, t_max):
        out_h = [0]
        i_prev_cur = 0
        i_h = gl_out[t]
        for i_cur in i_h:
            if i_cur[-3] > i_prev_cur:
                out_h.append(i_cur[-1])
                i_prev_cur = i_cur[-3]
        out_h.append(i_h[-1][-1] + 1)
        out_next_l.append(out_h)
    return out_next_l


def get_joint_prob(nu_list, nu_index_full, t_ind):
    ih = nu_index_full[t_ind]
    out_prob = []
    for i_full in ih:
        p_here = nu_list[0][1][i_full[0]]
        for t in range(1, t_ind + 1):
            ind = nu_list[t][i_full[2 * (t - 1)]][0].index(i_full[2 * t])
            p_here *= nu_list[t][i_full[2 * (t - 1)]][1][ind]
        out_prob.append(p_here)
    return out_prob


# FUNCTION THAT PRODUCES ADAPTED EMPIRICAL MEASURE USING KMEANS
# Basically implements the methodology of Backhoff et al "Estimating processes in adapted Wasserstein distance"
# However, instead of the fixed grid based approach therein, we use K-Means to cluster the points together, which is
# more flexible with respect to varying orders of magnitude etc
def empirical_k_means_measure(
    data, use_klist=0, klist=(), tol_decimals=6, use_weights=0, heuristic=0
):
    # data is [k, T_h] array
    # klist is list with T_h entries, each being an integer lower than k; number of barycenters for each time step
    (k, T_h) = data.shape
    if not use_klist:
        klist = (np.ones(T_h) * int(np.round(np.sqrt(k)))).astype(int)

    label_list = []
    support_list = []
    out_x = np.zeros([0, T_h])
    out_w = []

    # cluster points at each time point
    # print('Clustering...')
    if heuristic:
        for t in range(T_h):
            data_t = data[:, t]
            inds_sort_t = np.argsort(data_t)
            datas_t = data_t[inds_sort_t]
            n_av = int(np.round(k / klist[t]))
            lmax = int(np.floor(n_av * klist[t]))
            all_but_end = np.reshape(datas_t[:lmax], (-1, n_av))
            mean_all_but = np.mean(all_but_end, axis=1, keepdims=1)
            cx = mean_all_but
            mean_all_but = np.tile(mean_all_but, (1, n_av))
            mean_all_but = np.reshape(mean_all_but, (-1, 1))
            mean_rest = np.mean(datas_t[lmax:])
            if lmax < k:
                mean_vec = np.concatenate(
                    [np.squeeze(mean_all_but), np.array([mean_rest])]
                )
                cx = np.concatenate([cx, np.array([mean_rest])])
            else:
                mean_vec = np.squeeze(mean_all_but)
            lx = np.zeros(k, dtype=int)
            for i in range(k):
                for j in range(len(cx)):
                    if mean_vec[inds_sort_t[i]] == cx[j]:
                        lx[i] = j
                        continue
            label_list.append(lx)
            support_list.append(cx)

    else:
        for t in range(T_h):
            # print('t = ' + str(t))
            data_t = data[:, t : t + 1]
            kmx = KMeans(n_clusters=klist[t]).fit(data_t)
            cx = kmx.cluster_centers_
            cx = np.round(cx, decimals=tol_decimals)
            lx = kmx.labels_
            label_list.append(lx)
            support_list.append(cx)

    if (
        use_weights == 0
    ):  # weight all cluster centers equally? ... Convenient but theoretically flawed I think
        out = np.zeros([k, T_h])
        for t in range(T_h):
            out[:, t] = support_list[t][label_list[t]][:, 0]
        return out

    # build output measure
    for i in range(k):
        cur_path = np.zeros(T_h)
        for t in range(T_h):
            cur_path[t] = support_list[t][label_list[t][i]]

        # check whether the path already exists
        path_is_here = 0
        for j in range(len(out_w)):
            if np.all(out_x[j, :] == cur_path):
                out_w[j] += 1 / k
                path_is_here = 1
                break
        if not path_is_here:
            out_x = np.append(out_x, np.expand_dims(cur_path, axis=0), axis=0)
            out_w.append(1 / k)

    return out_x, out_w


# FUNCTIONS TO GENERATE VARIOUS SIMPLE PROBABILITY MEASURES
def brownian_motion_sample(t, s):
    # t is number of time steps
    # s is size of sample
    # returns brownian motion sample of size s x (t+1), where the zeroth time step is just constant zero
    incs = np.random.normal(size=[s, t])
    bms = np.cumsum(incs, axis=1)
    bms = np.concatenate([np.zeros([s, 1]), bms], axis=1)
    return bms


def binomial():
    # returns a measure which corresponds to a binomial model on a T time steps specified corresponding to its
    # natural graph, which is the temporal markovian graph.

    def mu(node, x_parents):
        if node == 0:
            return [[0], [1]]
        else:
            x = x_parents[
                0
            ]  # should only contain one element as the structure is Markovian
            return [[x - 1, x + 1], [0.5, 0.5]]

    def sup_mu(node_list):
        if len(node_list) == 0:
            out = np.array([])
            out = out.reshape(-1, 1)
            return out
        if len(node_list) == 1:
            t = node_list[
                0
            ]  # we only need to supply support for single element node-lists, as no two nodes share a child
            sup_l = [x for x in range(-t, t + 1, 2)]
        else:
            lol = []
            node_list_c = sorted(list(node_list))
            for t in node_list_c:
                lol.append([x for x in range(-t, t + 1, 2)])
            sup_l = list(product(*lol))
        bt = len(node_list)
        return np.reshape(np.array(sup_l), (-1, bt))

    return mu, sup_mu


def trinomial():
    # returns a measure which corresponds to a binomial model on a T time steps specified corresponding to its
    # natural graph, which is the temporal markovian graph.

    def mu(node, x_parents):
        if node == 0:
            return [[0], [1]]
        else:
            x = x_parents[
                0
            ]  # should only contain one element as the structure is Markovian
            return [[x - 1, x, x + 1], [1 / 3, 1 / 3, 1 / 3]]

    def sup_mu(node_list):
        if len(node_list) == 0:
            out = np.array([])
            out = out.reshape(-1, 1)
            return out
        if len(node_list) == 1:
            t = node_list[
                0
            ]  # we only need to supply support for single element node-lists, as no two nodes share a child
            sup_l = [x for x in range(-t, t + 1)]
        else:
            lol = []
            node_list_c = sorted(list(node_list))
            for t in node_list_c:
                lol.append([x for x in range(-t, t + 1)])
            sup_l = list(product(*lol))
        bt = len(node_list)
        return np.reshape(np.array(sup_l), (-1, bt))

    return mu, sup_mu


def rand_tree_pichler(T, num_branch=(2, 3, 2, 3, 4), init=10, udrange=10, discr=0):
    # Trying to implement the method in Chapter 5 from https://arxiv.org/pdf/2102.05413.pdf

    transitions = (
        {}
    )  # take as key a tuple of an (integer node and integer value) and returns a measure
    supports = {}  # takes as input a node and returns a set of integer support points
    for i in range(T + 1):
        supports[i] = set([])

    for t in range(T):
        if t == 0:
            if discr == 0:
                supports[0] = {init}
            else:
                supports[0] = {0}
            nbh = num_branch[0]
            probs = np.random.random_sample(nbh)
            probs = probs / np.sum(probs)
            if not discr:
                supps = np.random.randint(-udrange, udrange, size=[nbh, 1])
                supps_int = set(np.squeeze(10 + supps, axis=1))
                supports[1] |= supps_int
                transitions[(t + 1, init)] = [10 + supps, probs]
            else:
                supps = np.arange(0, nbh)
                supps = supps.reshape(-1, 1)
                supps_int = set(np.squeeze(supps, axis=1))
                supports[1] |= supps_int
                transitions[(t + 1, 0)] = [supps, probs]
        else:
            for x_int in supports[t]:
                nbh = num_branch[t]
                probs = np.random.random_sample(nbh)
                probs = probs / np.sum(probs)
                if not discr:
                    supps = np.random.randint(-udrange, udrange, size=[nbh, 1])
                    supps_int = set(np.squeeze(x_int + supps, axis=1))
                    supports[t + 1] |= supps_int
                    transitions[(t + 1, x_int)] = [x_int + supps, probs]
                else:
                    supps = np.arange(0, nbh)
                    supps = supps.reshape(-1, 1)
                    supps_int = set(np.squeeze(supps, axis=1))
                    supports[t + 1] |= supps_int
                    transitions[(t + 1, x_int)] = [supps, probs]

    if discr == 0:

        def mu(node, x_parents):
            if node == 0:
                return [np.reshape(np.array([10]), (-1, 1)), [1]]
            x = x_parents[
                0
            ]  # should only contain one element as the structure is Markovian
            x = int(x)
            return transitions[(node, x)]

    else:

        def mu(node, x_parents):
            if node == 0:
                return [np.reshape(np.array([0]), (-1, 1)), [1]]
            x = x_parents[
                0
            ]  # should only contain one element as the structure is Markovian
            x = int(x)
            return transitions[(node, x)]

    def sup_mu(node_list):
        if len(node_list) == 0:
            out = np.array([])
            out = out.reshape(-1, 1)
            return out
        return np.reshape(
            np.array(list(supports[node_list[0]])), (-1, 1)
        )  # we only supply support for single nodes

    print("Warning: You are using a measure where only one-step supports are specified")
    return mu, sup_mu


# CLASS FOR GRAPHS FOR REPRESENTING FUNCTIONS
# taken from https://www.geeksforgeeks.org/python-program-for-topological-sorting/
# Class to represent a graph
class Graph:
    def __init__(self, vertices):
        self.graph = defaultdict(list)  # dictionary containing adjacency List
        self.V = vertices  # No. of vertices
        self.parents = defaultdict(list)

    # function to add an edge to graph
    def addEdge(self, u, v):
        self.graph[u].append(v)
        self.parents[v].append(u)

        # A recursive function used by topologicalSort

    def topologicalSortUtil(self, v, visited, stack):

        # Mark the current node as visited.
        visited[v] = True

        # Recur for all the vertices adjacent to this vertex
        for i in self.graph[v]:
            if visited[i] == False:
                self.topologicalSortUtil(i, visited, stack)

                # Push current vertex to stack which stores result
        stack.insert(0, v)

        # The function to do Topological Sort. It uses recursive

    # topologicalSortUtil()
    def topologicalSort(self):
        # Mark all the vertices as not visited
        visited = [False] * self.V
        stack = []

        # Call the recursive helper function to store Topological
        # Sort starting from all vertices one by one
        for i in range(self.V):
            if visited[i] == False:
                self.topologicalSortUtil(i, visited, stack)

                # Print contents of stack
        return stack


# FUNCTION TO CHANGE REPRESENTATION OF MEASURES
def reduce_meas(marg, filt=0):
    # marg is a list with two entries, one n x d np array, and one n-list of weights
    # the goal of this function is to identify duplicates in the first entry and reduce the representation
    if filt == 1:
        return marg  # TODO: check how to reduce in this case as well ...
    if len(marg[0].shape) == 1:
        marg[0].reshape(-1, 1)
    uniques, inv = np.unique(marg[0], axis=0, return_inverse=1)
    w_new = np.zeros(len(uniques))
    for l in range(len(marg[0])):
        w_new[inv[l]] += marg[1][l]
    # for i in range(len(uniques)):
    #     for l in range(len(marg[0])):
    #         if np.all(uniques[i, :] == marg[0][l, :]):
    #             w_new[i] += marg[1][l]
    return [uniques, w_new]


def get_marginals_multi(
    mu, supp_mu, node_set, g, given_margs, all_out=0, index_mu=0, tol=10**-6, filt=0
):
    # function should get a joint distribution on the specified node set. I.e., we want to go from disintegration
    # representation of a measure towards specifying the joint distribution.
    # node_set is a tuple containing the nodes that we wish to calculate the joint marginal on

    # get relevant parents
    rel_par = []
    for nh in node_set:
        rel_par.extend(g.parents[nh].copy())
    rel_par = list(tuple(rel_par))
    rel_par.sort()
    rel_par = tuple(rel_par)
    rel_par_arr = np.array(rel_par)

    # for each node in node_set, get the indices of the respective parents in rel_par
    # rel_par_arr[indices[j]] will give the parents of node_set[j]
    indices = []
    for nh in node_set:
        ph = g.parents[nh].copy()
        index = np.zeros(len(ph), dtype=int)
        for ind0, i in enumerate(ph):
            for ind in range(len(rel_par)):
                if rel_par[ind] == i:
                    index[ind0] = ind
        indices.append(index)

    # get relevant marginal rel_marg of the form [ n x d array, list of weights]
    if len(rel_par) == 0:
        # for each node in node_set, get the parent values from xh
        x_lists = []
        w_lists = []
        if filt == 1:
            filt_list_a = []
            filt_list_b = []
        for ind, nh in enumerate(node_set):
            if index_mu == 0:
                marg_x_nh = mu(nh, [])
            else:
                marg_x_nh = mu(nh, 0)
            if len(np.array(marg_x_nh[0]).shape) == 2:
                x_lists.append(
                    marg_x_nh[0][:, 0]
                )  # marg_x_nh[0] is always of shape [n, 1]
                if filt == 1:
                    filt_list_a.append(marg_x_nh[2][:, 0])
                    filt_list_b.append(marg_x_nh[2][:, 1])
            else:
                x_lists.append(marg_x_nh[0])
                if filt == 1:
                    filt_list_a.append(marg_x_nh[2][:, 0])
                    filt_list_b.append(marg_x_nh[2][:, 1])

            w_lists.append(marg_x_nh[1])

        x_list_comb = list(product(*x_lists))
        if filt == 1:
            filt_a_comb = list(product(*filt_list_a))
            filt_b_comb = list(product(*filt_list_b))
            a_arr = np.array(filt_a_comb)
            b_arr = np.array(filt_b_comb)
            if len(a_arr.shape) == 1:
                a_arr = a_arr.reshape(-1, 1)
            if len(b_arr.shape) == 1:
                b_arr = b_arr.reshape(-1, 1)
            filt_tot = np.append(a_arr, b_arr, axis=1)
            w_list_comb = list(product(*w_lists))
            w_list_comb = [np.product(w_here) for w_here in w_list_comb]
            marg_add = [np.array(x_list_comb), w_list_comb, filt_tot]
        else:
            w_list_comb = list(product(*w_lists))
            w_list_comb = [np.product(w_here) for w_here in w_list_comb]
            marg_add = [np.array(x_list_comb), w_list_comb]
            marg_add = reduce_meas(marg_add, filt=filt)

        if all_out == 0:
            return marg_add
        else:
            node_set = list(node_set)
            node_set.sort()
            given_margs[tuple(node_set)] = marg_add
            return given_margs
    elif not rel_par in given_margs:
        if all_out == 0:
            rel_marg = get_marginals_multi(
                mu, supp_mu, rel_par, g, given_margs, index_mu=index_mu, filt=filt
            )
        else:
            all_margs = get_marginals_multi(
                mu,
                supp_mu,
                rel_par,
                g,
                given_margs,
                all_out=1,
                index_mu=index_mu,
                filt=filt,
            )
            rel_marg = all_margs[rel_par]
    else:
        rel_marg = given_margs[rel_par]
        all_margs = given_margs

    # calculate the joint marginal of node_set given the joint marginal of all the relevant parents:
    d = len(node_set)
    out_x = np.zeros([0, d])
    out_w = []
    if filt == 1:
        out_f = np.zeros([0, 2 * d])
    if index_mu == 1:
        if filt == 1:
            supp_rel_p = supp_mu(rel_par, filt=1)
            supp_rel = supp_rel_p[0]
            filt_rel = supp_rel_p[1]
        else:
            supp_rel = supp_mu(rel_par)
    for i in range(len(rel_marg[0])):
        xh = rel_marg[0][i, :]
        wh = rel_marg[1][i]
        if filt == 1:
            fh = rel_marg[2][
                i, :
            ]  # filtration of marginal; second dimension should be double that of xh
        if wh == 0:
            continue

        if index_mu == 1:
            ind_rel_here = -1
            for j in range(len(supp_rel)):
                if filt == 0:
                    if np.all(np.abs(xh - supp_rel[j, :]) < tol):
                        ind_rel_here = j
                        break
                else:
                    if np.all(np.abs(xh - supp_rel[j, :]) < tol) and np.all(
                        np.abs(fh - filt_rel[j, :]) < tol
                    ):
                        ind_rel_here = j
                        break
            if ind_rel_here == -1:
                print("ERROR: relevant support point not found...")

        # for each node in node_set, get the parent values from xh
        x_lists = []
        w_lists = []
        if filt == 1:
            filt_list_a = []
            filt_list_b = []
        for ind, nh in enumerate(node_set):
            if nh in rel_par:
                ind_nh = np.where(nh == rel_par_arr)[0][0]
                x_lists.append([xh[ind_nh]])
                w_lists.append([1])
                if filt == 1:
                    filt_list_a.append([fh[ind_nh]])
                    filt_list_b.append([fh[len(rel_par) + ind_nh]])
            else:
                if len(indices[ind]) > 0:
                    rel_x = xh[indices[ind]]
                    if filt == 1:
                        rel_filt_a = fh[indices[ind]]
                        rel_filt_b = fh[len(rel_par) + indices[ind]]

                else:
                    rel_x = []
                    if filt == 1:
                        rel_filt_a = []
                        rel_filt_b = []
                if index_mu == 0:
                    marg_x_nh = mu(nh, rel_x)
                else:
                    marg_x_nh = mu(nh, ind_rel_here)
                if len(np.array(marg_x_nh[0]).shape) == 2:
                    x_lists.append(
                        marg_x_nh[0][:, 0]
                    )  # marg_x_nh[0] is always of shape [n, 1]
                    if filt == 1:
                        filt_list_a.append(marg_x_nh[2][:, 0])
                        filt_list_b.append(marg_x_nh[2][:, 1])
                else:
                    x_lists.append(marg_x_nh[0])
                    if filt == 1:
                        filt_list_a.append(marg_x_nh[2][:, 0])
                        filt_list_b.append(marg_x_nh[2][:, 1])

                w_lists.append(marg_x_nh[1])

        if filt == 1:
            filt_a_comb = list(product(*filt_list_a))
            filt_b_comb = list(product(*filt_list_b))
            a_arr = np.array(filt_a_comb)
            b_arr = np.array(filt_b_comb)
            if len(a_arr.shape) == 1:
                a_arr = a_arr.reshape(-1, 1)
            if len(b_arr.shape) == 1:
                b_arr = b_arr.reshape(-1, 1)
            filt_tot = np.append(a_arr, b_arr, axis=1)
            x_list_comb = list(product(*x_lists))
            w_list_comb = list(product(*w_lists))
            w_list_comb = [np.product(w_here) * wh for w_here in w_list_comb]
            marg_add = [np.array(x_list_comb), w_list_comb, filt_tot]
            out_x = np.append(out_x, marg_add[0], axis=0)
            out_w.extend(marg_add[1])
            out_f = np.append(out_f, filt_tot, axis=0)
        else:
            x_list_comb = list(product(*x_lists))
            w_list_comb = list(product(*w_lists))
            w_list_comb = [np.product(w_here) * wh for w_here in w_list_comb]
            marg_add = [np.array(x_list_comb), w_list_comb]
            marg_add = reduce_meas(marg_add, filt=filt)
            out_x = np.append(out_x, marg_add[0], axis=0)
            out_w.extend(marg_add[1])

    if filt == 1:
        marg_out = [out_x, out_w, out_f]
    else:
        marg_out = [out_x, out_w]
    if all_out == 0:
        return reduce_meas(marg_out, filt=filt)
    else:
        node_set = list(node_set)
        node_set.sort()
        if tuple(node_set) not in all_margs:
            all_margs[tuple(node_set)] = reduce_meas(marg_out, filt=filt)
        for key in given_margs:
            if key not in all_margs:
                all_margs[key] = given_margs[key]
        return all_margs
