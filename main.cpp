#include <iostream>
#include <vector>
#include <Eigen/Dense>

Eigen::MatrixXd Lmatrix2paths(
    Eigen::MatrixXd L,
    int n_sample,
    bool normalize = false,
    int seed = 0,
    bool verbose = false
);

double Nested(Eigen::MatrixXd& X, Eigen::MatrixXd& Y);

int main() {

    Eigen::MatrixXd L(3, 3);
    L << 1, 0, 0,
         2, 4, 0,
         3, 2, 1;

    Eigen::MatrixXd M(3, 3);
    M << 1, 0, 0,
         2, 3, 0,
         3, 1, 2;

    int n_sample = 10000;
    Eigen::MatrixXd X = Lmatrix2paths(L, n_sample, true, 0, true);
    Eigen::MatrixXd Y = Lmatrix2paths(M, n_sample, true, 0, true);

    double res = Nested(X, Y);
}

