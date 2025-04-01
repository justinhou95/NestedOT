#include "header.h"

// Printer
void print_kernel_x(std::vector<ConditionalDistribution>& kernel_x);

// Generate Data Structure
Eigen::MatrixXd Lmatrix2paths(
    Eigen::MatrixXd L,
    int n_sample,
    bool normalize = false,
    int seed = 0,
    bool verbose = false
);

Eigen::MatrixXd path2adaptedpath(const Eigen::MatrixXd& samples, double delta_n);

void v_set_add(const Eigen::MatrixXd& mat, std::set<double>& unique_set);

Eigen::MatrixXi quantize_path(Eigen::MatrixXd& adaptedX, std::map<double, int>& v2q);

Eigen::MatrixXi sort_qpath(const Eigen::MatrixXi& path);

std::vector<std::map<std::vector<int>, std::map<int, int>>> qpath2mu_x(Eigen::MatrixXi& qpath);

std::vector<ConditionalDistribution> mu_x2kernel_x(std::vector<std::map<std::vector<int>, std::map<int, int>>>& mu_x);

// Solver
std::vector<std::vector<float>> SquareCost(std::vector<int>& vx, std::vector<int>& vy, std::vector<std::vector<float>>& cost_matrix);

void AddDppValue(std::vector<std::vector<float>>& cost, std::vector<std::vector<float>>& Vtplus, int& i0, int& j0);

float SolveOT(std::vector<float>& wx, std::vector<float>& wy, std::vector<std::vector<float>>& cost);
