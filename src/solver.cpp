#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <set>
#include <map>
#include "header_dist.h"


Eigen::MatrixXd path2adaptedpath(const Eigen::MatrixXd& samples, double grid_size);

void v_set_add(const Eigen::MatrixXd& mat, std::set<double>& unique_set);

Eigen::MatrixXi quantize_path(Eigen::MatrixXd& adaptedX, std::map<double, int>& v2q);

Eigen::MatrixXi sort_qpath(const Eigen::MatrixXi& path);

std::vector<std::map<std::vector<int>, std::map<int, int>>> qpath2mu_x(Eigen::MatrixXi& qpath, const bool& markovian);

std::vector<ConditionalDistribution> mu_x2kernel_x(std::vector<std::map<std::vector<int>, std::map<int, int>>>& mu_x);


int EMD_wrap(int n1, int n2, double *X, double *Y, double *D, double *G,
    double* alpha, double* beta, double *cost, uint64_t maxIter);

std::vector<std::vector<double>> SquareCost(std::vector<int>& vx, std::vector<int>& vy, std::vector<std::vector<double>>& cost_matrix){
    int m = vx.size();
    int n = vy.size();
    std::vector<std::vector<double>> cost(m, std::vector<double>(n, 0.0f));
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            cost[i][j] = cost_matrix[vx[i]][vy[j]];
        }
    }
    return cost;
}

void AddDppValue(std::vector<std::vector<double>>& cost, std::vector<std::vector<double>>& Vtplus, int& i0, int& j0){
    int m = cost.size();
    int n = cost[0].size();
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            cost[i][j] += Vtplus[i0+i][j0+j];
        }
    }
}

void AddDppValueMarkovian(
    std::vector<std::vector<double>>& cost, 
    std::vector<std::vector<double>>& Vtplus, 
    std::vector<int>& x_next_idx,
    std::vector<int>& y_next_idx){
    int m = cost.size();
    int n = cost[0].size();
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            cost[i][j] += Vtplus[x_next_idx[i]][y_next_idx[j]];
        }
    }
}


double SolveOT(std::vector<double>& wx, std::vector<double>& wy, std::vector<std::vector<double>>& cost){
    int n1 = wx.size(); 
    int n2 = wy.size();
    double c = 0.0;

    if (n1 == 1 || n2 == 1){
        for (int i = 0; i < n1; i++){
            for (int j = 0; j < n2; j++){
                c += wx[i] * cost[i][j] * wy[j];
            }
        }
    } else {
        double* X = wx.data(); // First marginal histogram 
        double* Y = wy.data(); // Second marginal histogram
        double* C = new double[n1 * n2]; // Cost matrix
        for (int i = 0; i < n1; ++i)
            for (int j = 0; j < n2; ++j)
                C[i * n2 + j] = cost[i][j];
        double* G = new double[n1 * n2]; // Coupling
        double* alpha = new double[n1]; // First dual potential function  
        double* beta = new double[n2]; // Second dual potential function          
        uint64_t maxIter = 100000;
        
        int result = EMD_wrap(n1, n2, X, Y, C, G, alpha, beta, &c, maxIter);

        if (result != 1){
            std::cout << "OT is not solved optimally" << std::endl;
        } 

        delete[] C;
        delete[] G;
        delete[] alpha;
        delete[] beta;
    }

    return c;

}




// Parameters:
//   X, Y         : Input paths (each row is a time step, columns are samples)
//   grid_size    : Grid size used for adapting/quantizing the paths.
//   markovian    : Switch between markovian (true) and full history (false) processing.
//   num_threads  : Number of threads to use (if <= 0, maximum available threads are used).
//   power        : Exponent for the cost function (only power 1 and 2 are optimized here).

double Nested(Eigen::MatrixXd& X,
    Eigen::MatrixXd& Y,
    double grid_size,
    const bool& markovian,
    int num_threads,
    const int power,
    const bool verbose)
{

    // Determine number of threads.
    if (num_threads <= 0) {
        num_threads = omp_get_max_threads();
    }

    int T = X.rows()-1;
    int n_sample = X.cols();
    Eigen::MatrixXd adaptedX = path2adaptedpath(X, grid_size);
    Eigen::MatrixXd adaptedY = path2adaptedpath(Y, grid_size);

    std::set<double> v_set;
    v_set_add(adaptedX, v_set);
    v_set_add(adaptedY, v_set);

    std::map<double, int> v2q; // Map value to quantization e.g. v2q[3.5] = 123
    std::vector<double> q2v;  // Map quantization to value e.g.  q2v[123] = 3.5
    int pos = 0;
    for (double v : v_set) {
        v2q[v] = pos;
        q2v.push_back(v);
        pos += 1;
    }

    Eigen::MatrixXi qX = sort_qpath(quantize_path(adaptedX, v2q).transpose());
    Eigen::MatrixXi qY = sort_qpath(quantize_path(adaptedY, v2q).transpose());


    // std::cout << qX << std::endl;

    std::vector<std::map<std::vector<int>, std::map<int, int>>> mu_x = qpath2mu_x(qX, markovian);
    std::vector<std::map<std::vector<int>, std::map<int, int>>> nu_y = qpath2mu_x(qY, markovian);

    // for (int t = 0; t < T; t++){
    //     print_mu_x(mu_x[t]);
    // }   

    std::vector<ConditionalDistribution> kernel_x = mu_x2kernel_x(mu_x);
    std::vector<ConditionalDistribution> kernel_y = mu_x2kernel_x(nu_y);

    // print_kernel_x(kernel_x);
    if (verbose){std::cout << "Start computing with " << num_threads << " threads" << std::endl;}

    std::vector<std::vector<std::vector<double>>> V(T);
    for(int t=0; t<T; t++){
        V[t] = std::vector<std::vector<double>>(kernel_x[t].nc, std::vector<double>(kernel_y[t].nc, 0.0f));
    }


    // Preselect the cost function.
    // For performance we handle only power 1 and 2 here (the inner loop will only access the cost matrix).
    std::function<double(double)> cost_func;
    if (power == 1) {
        cost_func = [](double diff) { return std::abs(diff); };
    } else if (power == 2) {
        cost_func = [](double diff) { return diff * diff; };
    } else {
        // Fallback (this branch will be used infrequently).
        cost_func = [power](double diff) { return std::pow(std::abs(diff), power); };
    }

    // Precompute the cost matrix (base cost between quantized values).
    std::vector<std::vector<double>> cost_matrix(q2v.size(), std::vector<double>(q2v.size(), 0.0));
    for (int i = 0; i < (int)q2v.size(); i++) {
        for (int j = 0; j < (int)q2v.size(); j++) {
            double diff = q2v[i] - q2v[j];
            cost_matrix[i][j] = cost_func(diff);
        }
    }


    auto start = std::chrono::steady_clock::now();

    for (int t = T - 1; t >= 0; t--){
        if (verbose){
        std::cout << "Timestep " << t << std::endl;
        std::cout << "Computing " <<  kernel_x[t].nc * kernel_y[t].nc << " OTs ......." << std::endl;
        }
        #pragma omp parallel for num_threads(num_threads) if(kernel_x[t].nc > 100)
        for (int ix = 0; ix < kernel_x[t].nc; ix++){
            for (int iy = 0; iy < kernel_y[t].nc; iy++){
                // Reference with shorter names
                std::vector<int>& vx = kernel_x[t].dists[ix].values;
                std::vector<int>& vy = kernel_y[t].dists[iy].values;
                std::vector<double>& wx = kernel_x[t].dists[ix].weights;
                std::vector<double>& wy = kernel_y[t].dists[iy].weights;
                
                std::vector<std::vector<double>> cost = SquareCost(vx, vy, cost_matrix);
                if (t < T - 1){
                    if (markovian){
                        std::vector<int>& x_next_idx = kernel_x[t].next_idx[ix];
                        std::vector<int>& y_next_idx = kernel_y[t].next_idx[iy];
                        AddDppValueMarkovian(cost, V[t + 1], x_next_idx, y_next_idx);
                    } else {
                        int& i0 = kernel_x[t].nv_cums[ix];
                        int& j0 = kernel_y[t].nv_cums[iy];
                        AddDppValue(cost, V[t + 1], i0, j0);
                    }
                }
                V[t][ix][iy] = SolveOT(wx, wy, cost);
            }
        }
    }

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << std::chrono::duration<double, std::milli>(diff).count()/1000. << " seconds" << std::endl;

    double nested_ot_value = V[0][0][0];
    if (verbose){
        std::cout << "Nested OT value: " << nested_ot_value << std::endl;
        std::cout << "Finish" << std::endl;
    }
    

    return nested_ot_value;
}