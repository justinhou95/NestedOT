import matplotlib.pyplot as plt
import numpy as np


def nested(X, Y, delta_n, markovian):
    try:
        import _wrapper

        return _wrapper.nested(X, Y, delta_n, markovian)
    except:
        pass
        # python_nested(X, Y, delta_n, markovian)


def nestedmarkovian(X, Y, grid_size):
    import _wrapper

    return _wrapper.nestedmarkovian(X, Y, grid_size)
