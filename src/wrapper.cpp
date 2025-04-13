#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <Eigen/Dense>
#include <vector>

namespace py = pybind11;

double Nested(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, double& delta_n, const bool& markovian);

double NestedBindToPython(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, double& delta_n, const bool& markovian){
    py::scoped_ostream_redirect stream(
        std::cout,                               // std::ostream&
        py::module_::import("sys").attr("stdout") // Python stdout
    );
    return Nested(X, Y, delta_n, markovian);
}


PYBIND11_MODULE(_wrapper, m) {
    m.doc() = R"pbdoc(pnot)pbdoc";

    m.def("nested_ot_solver", &NestedBindToPython, R"pbdoc(
        Nested OT
    )pbdoc");

}