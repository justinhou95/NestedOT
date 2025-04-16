from warnings import WarningMessage


def nested_ot(X, Y, grid_size, markovian, parallel=True, num_threads=8, power=2):
    try:
        import _wrapper

        if not parallel:
            WarningMessage("Using C++ solver (always parallel)")

        return _wrapper.nested_ot_solver(
            X, Y, grid_size, markovian, num_threads, power
        )  # Always run in parallel as the goal is speed with C++ solver
    except:
        WarningMessage(
            "Can not find the C++ solver, using the Python one (much slower)"
        )
        from .solver import nested_ot_solver_py

        print("Using Python solver")

        return nested_ot_solver_py(
            X, Y, grid_size, markovian, parallel, num_threads, power
        )
