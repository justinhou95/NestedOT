#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <set>
#include <map>
#include "header.h"
#include "header_main.h"


int EMD_wrap(int n1, int n2, double *X, double *Y, double *D, double *G,
    double* alpha, double* beta, double *cost, uint64_t maxIter);

int main() {
    // Test the code LP solver
    int n1 = 3, n2 = 3;
    double X1[] = {0.3, 0.4, 0.3};  // Distribution 1 (supply)
    double Y1[] = {0.2, 0.5, 0.3};  // Distribution 2 (demand)

    // Distance matrix D[n1][n2] in row-major order
    double D[] = {
        0.0, 1.0, 2.0,
        1.0, 0.0, 1.0,
        2.0, 1.0, 0.0
    };

    double G[9] = {0};            // Flow matrix (output)
    double alpha[3] = {0};          // Dual potentials (output)
    double beta[3] = {0};           // Dual potentials (output)
    double cost = 0.0;

    uint64_t maxIter = 100000;

    int result = EMD_wrap(n1, n2, X1, Y1, D, G, alpha, beta, &cost, maxIter);

    if (result == 1) {
        std::cout << "EMD cost: " << cost << std::endl;
        std::cout << "Flow matrix G:" << std::endl;
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                std::cout << G[i * n2 + j] << " ";
            }
            std::cout << std::endl;
        }
    } else {
        std::cerr << "EMD computation failed with codeeee: " << result << std::endl;
    }

    // Compute numerical adapted Wasserstein distance

    constexpr int T = 3;
    Eigen::MatrixXd L(3, 3);
    L << 1, 0, 0,
         2, 4, 0,
         3, 2, 1;

    Eigen::MatrixXd M(3, 3);
    M << 1, 0, 0,
         2, 3, 0,
         3, 1, 2;

    int n_sample = 10;
    Eigen::MatrixXd X = Lmatrix2paths(L, n_sample, true, 0, true);
    Eigen::MatrixXd Y = Lmatrix2paths(M, n_sample, true, 0, true);

    double delta_n = 0.1;
    Eigen::MatrixXd adaptedX = path2adaptedpath(X, delta_n);
    Eigen::MatrixXd adaptedY = path2adaptedpath(Y, delta_n);

    std::set<double> v_set;
    v_set_add(adaptedX, v_set);
    v_set_add(adaptedY, v_set);

    std::map<double, int> v2q; 
    std::vector<double> q2v; 

    int pos = 0;
    for (double v : v_set) {
        v2q[v] = pos;
        q2v.push_back(v);
        pos += 1;
    }

    Eigen::MatrixXi qX = sort_qpath(quantize_path(adaptedX, v2q).transpose());
    Eigen::MatrixXi qY = sort_qpath(quantize_path(adaptedY, v2q).transpose());

    std::vector<std::map<std::vector<int>, std::map<int, int>>> mu_x = qpath2mu_x(qX);
    std::vector<std::map<std::vector<int>, std::map<int, int>>> nu_y = qpath2mu_x(qX);

    std::vector<ConditionalDistribution> kernel_x = mu_x2kernel_x(mu_x);
    std::vector<ConditionalDistribution> kernel_y = mu_x2kernel_x(nu_y);

    // print_kernel_x(kernel_x);

    std::vector<std::vector<std::vector<double>>> V(T);
    for(int t=0; t<T; t++){
        V[t] = std::vector<std::vector<double>>(kernel_x[t].nc, std::vector<double>(kernel_y[t].nc, 0.0f));
    }

    std::vector<std::vector<double>> cost_matrix(q2v.size(), std::vector<double>(q2v.size(), 0.0f));
    for (int i = 0; i < q2v.size(); i++){
        for (int j = 0; j < q2v.size(); j++){
            cost_matrix[i][j] = (q2v[i] - q2v[j]) * (q2v[i] - q2v[j]);
        }
    }

    for (int t = T - 1; t >= 0; t--){
        std::cout << "Timestep " << t << std::endl;
        for (int ix = 0; ix < kernel_x[t].nc; ix++){
            for (int iy = 0; iy < kernel_y[t].nc; iy++){
                // Reference with shorter names
                std::vector<int>& vx = kernel_x[t].dists[ix].values;
                std::vector<int>& vy = kernel_y[t].dists[iy].values;
                std::vector<double>& wx = kernel_x[t].dists[ix].weights;
                std::vector<double>& wy = kernel_y[t].dists[iy].weights;
                int& i0 = kernel_x[t].nv_cums[ix];
                int& j0 = kernel_y[t].nv_cums[iy];

                std::vector<std::vector<double>> cost = SquareCost(vx, vy, cost_matrix);
                if(t < T-1) AddDppValue(cost, V[t+1], i0, j0);
                V[t][ix][iy] = SolveOT(wx, wy, cost);
            }
        }
    }

    double AW2 = V[0][0][0];
    std::cout << "AW_2^2: " << AW2 << std::endl;
    std::cout << "Finish" << std::endl;

    return 0;
}

