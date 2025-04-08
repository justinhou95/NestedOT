#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>


double Nested(Eigen::MatrixXd& X, Eigen::MatrixXd& Y);

namespace py = pybind11;


// OpenMP test
int sum_thread_ids(){
    int sum = 0;
    #pragma omp parallel for
    for (int i = 0; i<10; i++)
    {
        sleep(1);
    }
    return 0;
}

PYBIND11_MODULE(_wrapper, m) {
    m.doc() = R"pbdoc(Some Documentation about pnot)pbdoc";

    m.def("nested", &Nested, R"pbdoc(
        Some other explanation about the add function.
    )pbdoc");

    m.def("sum_thread_ids", &sum_thread_ids, "Adds the id of threads");
}