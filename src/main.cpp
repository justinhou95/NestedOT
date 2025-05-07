#include <iostream>
#include <vector>
#include <Eigen/Dense>

Eigen::MatrixXd Lmatrix2paths(
    Eigen::MatrixXd L,
    int n_sample,
    int seed = 0
);

double Nested(Eigen::MatrixXd& X,
    Eigen::MatrixXd& Y,
    double grid_size,
    const bool& markovian,
    int num_threads,
    const int power,
    const bool verbose);


double Nested2(Eigen::MatrixXd& X,
    Eigen::MatrixXd& Y,
    double grid_size,
    const bool& markovian,
    int num_threads,
    const int power,
    const bool verbose);

int main() {

    Eigen::MatrixXd L(3, 3);
    L << 1, 0, 0,
         2, 4, 0,
         3, 2, 1;

    Eigen::MatrixXd M(3, 3);
    M << 1, 0, 0,
         2, 3, 0,
         3, 1, 2;

    int n_sample = 80000;
    Eigen::MatrixXd X = Lmatrix2paths(L, n_sample, true);
    Eigen::MatrixXd Y = Lmatrix2paths(M, n_sample, true);
    const bool markovian = true;
    double delta_n = 0.01;


    double res = Nested(X, Y, delta_n, markovian, 8, 2, true);

    double res2 = Nested2(X, Y, delta_n, markovian, 8, 2, true);

    // bool markovians[2] = {true, false};
    // int num_threads[3] = {0, 1, 16};
    // int ps[3] = {1, 2, 3};

    // for (int p: ps){
    //     for (bool markovian: markovians){    
    //         for (int num_thread: num_threads){
    //             std::cout << "NumtThread: " << num_thread  << "  ploss: " << p << "  Markovian: " << markovian << std::endl;
    //             double res = Nested(X, Y, delta_n, markovian, num_thread, p, true);
    //             std::cout << "========================================================================" << std::endl;
    //         }
    //     }
    // }

}
