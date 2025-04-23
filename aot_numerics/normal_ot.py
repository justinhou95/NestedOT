import numpy as np
import ot

# This file contains different functions to solve a normal OT problem with two marginals using different methods.


def gurobi_2d(
    margs, f, p_dist=2, radial_cost=0, f_id=0, minmax="min", r_opti=0, outputflag=1
):
    """
    :param margs: list with 2 entries, each entry being a discrete probability measure
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
    n1 = len(xl_1)
    n2 = len(xl_2)

    if len(xl_1.shape) == 1:
        xl_1 = xl_1.reshape(-1, 1)
    if len(xl_2.shape) == 1:
        xl_2 = xl_2.reshape(-1, 1)

    # build cost matrix:
    if radial_cost == 0:
        cost_mat = np.zeros([n1, n2])
        for i in range(n1):
            for j in range(n2):
                cost_mat[i, j] = f(xl_1[i], xl_2[j])
    else:
        cost_mat = np.linalg.norm(
            xl_1[:, None, :] - xl_2[None, :, :], axis=-1, ord=p_dist
        )
        if f_id == 0:
            cost_mat = f(cost_mat)

    # initialize model
    m = Model("Primal")
    if outputflag == 0:
        m.setParam("OutputFlag", 0)
    pi_var = m.addVars(n1, n2, lb=0, ub=1, name="pi_var")

    # add marginal constraints
    m.addConstrs((pi_var.sum(i, "*") == pl_1[i] for i in range(n1)), name="first_marg")
    m.addConstrs((pi_var.sum("*", i) == pl_2[i] for i in range(n2)), name="second_marg")

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


def solve_unbalanced(
    margs, f, minmax="min", r_opti=0, outputflag=1, epsilon=0.5, alpha=10
):
    """
    should converge to normal ot for alpha --> infty and epsilon --> 0.
    In practice, epsilon lower than 0.01 may cause problems
    And on the other hand, larger values of epsilon and lower values of alpha converge a lot faster
    Notably, the penalization by alpha allows for more general couplings
    (that do not have to satisfy marginal constraints)
    On the other hand, regularization by epsilon simply restricts the couplings (by requiring smoothness)
    Hence for low values of alpha, optimal value is usually above the true OT
    And for high value sof epsilon, optimal value is usually below the true OT
    """
    m1 = margs[0]
    m2 = margs[1]
    xl_1 = np.array(m1[0])
    xl_2 = np.array(m2[0])
    pl_1 = m1[1]
    pl_2 = m2[1]
    n1 = len(xl_1)
    n2 = len(xl_2)

    if len(xl_1.shape) == 1:
        xl_1 = xl_1.reshape(-1, 1)
    if len(xl_2.shape) == 1:
        xl_2 = xl_2.reshape(-1, 1)

    cost_mat = np.zeros([n1, n2])
    for i in range(n1):
        for j in range(n2):
            cost_mat[i, j] = f(xl_1[i], xl_2[j])
    if minmax == "max":
        cost_mat *= -1
    Gs = ot.unbalanced.sinkhorn_unbalanced(
        pl_1, pl_2, cost_mat, epsilon, alpha, verbose=outputflag
    )

    if r_opti == 0:
        return np.sum(Gs * cost_mat)
    else:
        return np.sum(Gs * cost_mat), Gs


def solve_pot(margs, f, minmax="min", r_opti=0, outputflag=1):
    m1 = margs[0]
    m2 = margs[1]
    xl_1 = np.array(m1[0])
    xl_2 = np.array(m2[0])
    pl_1 = m1[1]
    pl_2 = m2[1]
    n1 = len(xl_1)
    n2 = len(xl_2)

    if len(xl_1.shape) == 1:
        xl_1 = xl_1.reshape(-1, 1)
    if len(xl_2.shape) == 1:
        xl_2 = xl_2.reshape(-1, 1)

    cost_mat = np.zeros([n1, n2])
    for i in range(n1):
        for j in range(n2):
            cost_mat[i, j] = f(xl_1[i], xl_2[j])
    if minmax == "max":
        cost_mat *= -1
    Gs = ot.emd(pl_1, pl_2, cost_mat)

    if r_opti == 0:
        return np.sum(Gs * cost_mat)
    else:
        return np.sum(Gs * cost_mat), Gs


def solve_sinkhorn(margs, f, minmax="min", r_opti=0, outputflag=1, epsilon=0.01):
    # should converge to normal ot for epsilon --> 0.
    # In practice, epsilon lower than around 0.01 may cause problems
    m1 = margs[0]
    m2 = margs[1]
    xl_1 = np.array(m1[0])
    xl_2 = np.array(m2[0])
    pl_1 = m1[1]
    pl_2 = m2[1]
    n1 = len(xl_1)
    n2 = len(xl_2)

    if len(xl_1.shape) == 1:
        xl_1 = xl_1.reshape(-1, 1)
    if len(xl_2.shape) == 1:
        xl_2 = xl_2.reshape(-1, 1)

    cost_mat = np.zeros([n1, n2])
    for i in range(n1):
        for j in range(n2):
            cost_mat[i, j] = f(xl_1[i], xl_2[j])
    if minmax == "max":
        cost_mat *= -1
    Gs = ot.sinkhorn(pl_1, pl_2, cost_mat, epsilon)

    if r_opti == 0:
        return np.sum(Gs * cost_mat)
    else:
        return np.sum(Gs * cost_mat), Gs
