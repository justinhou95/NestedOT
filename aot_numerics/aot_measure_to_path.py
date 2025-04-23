import numpy as np


# Function to transfer from "aot_numerics" framework to mine
def extract_sample_paths(mu, T, init):
    """
    Extract full sample paths and their weights from a disintegrated measure.

    Parameters:
      mu: A function representing the disintegrated measure. It accepts (node, x_parents)
          and returns a tuple [support, weights], where support is a 2D array (n x 1) of possible
          values at the given node and weights is a list (or array) of probabilities.
      T: Total number of time steps (i.e., nodes are 0, 1, ..., T). Note that T is the final node.
      init: The initial state (used at node 0).

    Returns:
      paths: A NumPy array of shape (N_paths, T+1) where each row is a sample path.
      weights: A NumPy array of shape (N_paths,) with the corresponding probability weights.
    """
    # Start at node 0 with the initial state.
    paths = [[init]]
    weights = [1.0]

    # For nodes 1 through T, extend each path using the measure mu.
    for t in range(1, T + 1):
        new_paths = []
        new_weights = []
        # For each path constructed so far:
        for path, w in zip(paths, weights):
            # For a Markov measure, we pass the last state as the only parent.
            parent = path[-1]
            # Get the disintegrated measure at node t.
            support, supp_weights = mu(t, [parent])
            # Each row in support is assumed to be an array of shape (1,); extract its value.
            for s, p in zip(support, supp_weights):
                # s is a 1D array (of length 1); extract its scalar value.
                new_val = s[0]
                new_paths.append(path + [new_val])
                new_weights.append(w * p)
        paths = new_paths
        weights = new_weights

    return np.array(paths), np.array(weights)
