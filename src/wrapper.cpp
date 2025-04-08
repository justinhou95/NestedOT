#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>


double Nested(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, double delta_n);

namespace py = pybind11;

PYBIND11_MODULE(_wrapper, m) {
    m.doc() = R"pbdoc(pnot)pbdoc";

    m.def("nested", &Nested, R"pbdoc(
        Nested OT
    )pbdoc");

}