from warnings import WarningMessage


def nested_ot(X, Y, delta_n, markovian, parallel=True):
    try:
        import _wrapper

        return _wrapper.nested_ot_solver(X, Y, delta_n, markovian)
    except:
        WarningMessage(
            "Can not find the C++ solver, using the Python one (much slower)"
        )
        from .py_solver import nested_ot_solver_py

        return nested_ot_solver_py(X, Y, delta_n, markovian, parallel)
