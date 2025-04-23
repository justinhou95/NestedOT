import numpy as np
from time import time
from sklearn.cluster import KMeans
from aot_numerics.normal_ot import *
from aot_numerics.measure import *
from aot_numerics.mainfunctions import *

"""
This script implements is to reproduce result in:

"Fitted Value Iteration Methods for Bicausal Optimal Transport"
by Erhan Bayraktar and Bingyan Han (2023).

Reference: https://arxiv.org/abs/2306.12658

The original code is sourced from: https://github.com/hanbingyan/FVIOT

The goal is to set up measures to reproduce the numerical experiments conducted in the paper.
"""


def rand_tree_binom(T, init, vol, N_leaf=2, in_size=100):
    # Trying to implement the method in Chapter 5 from https://arxiv.org/pdf/2102.05413.pdf

    transitions = (
        {}
    )  # take as key a tuple of an (integer node and integer value) and returns a measure
    supports = {}  # takes as input a node and returns a set of integer support points
    for i in range(T + 1):
        supports[i] = set([])

    for t in range(T):
        if t == 0:
            supports[0] = {init}
            rand_supps = np.random.randn(N_leaf * in_size, 1) * vol
            kmeans = KMeans(n_clusters=N_leaf, n_init=10).fit(rand_supps)
            _, probs = np.unique(kmeans.labels_, return_counts=True)
            probs = probs / np.sum(probs)
            supps = kmeans.cluster_centers_
            supps_int = set(np.squeeze(init + supps, axis=1))
            supports[1] |= supps_int
            transitions[(t + 1, init)] = [init + supps, probs]
        else:
            for x_int in supports[t]:
                rand_supps = np.random.randn(N_leaf * in_size, 1) * vol
                kmeans = KMeans(n_clusters=N_leaf, n_init=10).fit(rand_supps)
                _, probs = np.unique(kmeans.labels_, return_counts=True)
                # probs = np.ones(N_leaf)
                probs = probs / np.sum(probs)
                supps = kmeans.cluster_centers_
                supps_int = set(np.squeeze(x_int + supps, axis=1))
                # supps_int = set(x_int + supps)
                supports[t + 1] |= supps_int
                transitions[(t + 1, x_int)] = [x_int + supps, probs]

    def mu(node, x_parents):
        if node == 0:
            return [np.reshape(np.array([init]), (-1, 1)), [1]]
        x = x_parents[
            0
        ]  # should only contain one element as the structure is Markovian
        # x = int(x)
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


def comb_tree(data, T, init, klist):

    transitions = (
        {}
    )  # take as key a tuple of an (integer node and integer value) and returns a measure
    supports = {}  # takes as input a node and returns a set of integer support points
    for i in range(T + 1):
        supports[i] = set([])

    label_list = []
    support_list = []

    for t in range(T + 1):
        # print('t = ' + str(t))
        data_t = data[:, t : t + 1]
        kmx = KMeans(n_clusters=klist[t], n_init=10).fit(data_t)
        cx = kmx.cluster_centers_
        cx = np.round(cx, decimals=6)
        lx = kmx.labels_
        label_list.append(lx)
        support_list.append(cx)

    supports[0] = {init}
    for t in range(T):
        # if t == 0:
        #     supports[0] = {init}
        #     supps = support_list[1]
        #     _, probs = np.unique(label_list[0], return_counts=True)
        #     probs = probs/np.sum(probs)
        #     supps_int = set(np.squeeze(supps, axis=1))
        #     supports[1] |= supps_int
        #     transitions[(t+1, init)] = [supps, probs]
        # else:
        uniq_label = np.unique(label_list[t])
        for label_int in uniq_label:
            x_int = support_list[t][label_int][0]
            supps_labels, probs = np.unique(
                label_list[t + 1][label_list[t] == label_int], return_counts=True
            )
            probs = probs / np.sum(probs)
            supps = support_list[t + 1][supps_labels]
            # supps = set(np.squeeze(supps, axis=1))
            transitions[(t + 1, x_int)] = [supps, probs]
        supports[t + 1] = set(np.squeeze(support_list[t + 1], axis=1))

    def mu(node, x_parents):
        if node == 0:
            return [np.reshape(np.array([init]), (-1, 1)), [1]]
        x = x_parents[
            0
        ]  # should only contain one element as the structure is Markovian
        # x = int(x)
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
