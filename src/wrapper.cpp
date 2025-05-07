#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <Eigen/Dense>
#include <vector>

namespace py = pybind11;

// Forward-declare our Nested function.
double Nested(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, double grid_size, const bool& markovian, int num_threads, const int power, const bool verbose);

// A simple wrapper function to redirect C++ std::cout to Python's stdout.
double NestedPython(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, double grid_size, const bool& markovian, int num_threads = -1, const int power = 2, const bool verbose = false) {
    py::scoped_ostream_redirect stream(
        std::cout,                               // redirect C++ ostream
        py::module_::import("sys").attr("stdout") // to Python stdout
    );
    return Nested(X, Y, grid_size, markovian, num_threads, power, verbose);
}

PYBIND11_MODULE(_wrapper, m) {
    m.doc() = R"pbdoc(
        pnot OT Module
    )pbdoc";

    m.def("nested_ot_solver", &NestedPython, R"pbdoc(
        Nested optimal transport solver.
        Parameters:
          X, Y        : Input paths (Eigen::MatrixXd)
          grid_size   : Grid size (double)
          markovian   : Whether to use markovian processing (bool)
          num_threads : Number of threads (int; if <=0, max available threads are used)
          power       : Exponent for cost function (int, typically 1 or 2)
          verbose     : Whether print verbose information
        Returns:
          Adapted Wasserstein squared distance (double)
    )pbdoc");
}