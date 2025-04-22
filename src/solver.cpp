#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <chrono>
#include <set>
#include <omp.h> 
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



// Memory buffers for the EMD solver
static thread_local std::vector<double> D_buf, G_buf, α_buf, β_buf;

// flattened OT‐wrapper for Markovian case:
double SolveOT_flat(
    const std::vector<double>& wx,
    const std::vector<double>& wy,
    const std::vector<std::vector<double>>& cost_matrix,
    const std::vector<int>& vx,
    const std::vector<int>& vy,
    const std::vector<int>* x_next_idx,   // nullptr if t==T-1
    const std::vector<int>* y_next_idx,   // nullptr if t==T-1
    const std::vector<std::vector<double>>& Vtplus, // empty if t==T-1
    uint64_t maxIter = 100000
) {
    int n1 = wx.size(), n2 = wy.size();
    int sz = n1 * n2;

    // resize our buffers once per thread
    D_buf.resize(sz);
    G_buf.resize(sz);
    α_buf.resize(n1);
    β_buf.resize(n2);

    // fill D_buf as a single flat loop
    for (int i = 0; i < n1; ++i) {
        const int vi = vx[i];
        if (x_next_idx) {
            int xi = (*x_next_idx)[i];
            for (int j = 0; j < n2; ++j) {
                int vj = vy[j];
                int yj = (*y_next_idx)[j];
                D_buf[i*n2 + j]
                  = cost_matrix[vi][vj]
                  + Vtplus[xi][yj];
            }
        } else {
            for (int j = 0; j < n2; ++j) {
                int vj = vy[j];
                D_buf[i*n2 + j] = cost_matrix[vi][vj];
            }
        }
    }

    // trivial 1×N or N×1 case
    if (n1 == 1 || n2 == 1) {
        double c = 0.0;
        for (int i = 0; i < n1; ++i)
            for (int j = 0; j < n2; ++j)
                c += wx[i] * D_buf[i*n2 + j] * wy[j];
        return c;
    }

    // otherwise call EMD backend
    double cost_out = 0.0;
    int res = EMD_wrap(
      n1, n2,
      const_cast<double*>(wx.data()),
      const_cast<double*>(wy.data()),
      D_buf.data(), G_buf.data(),
      α_buf.data(), β_buf.data(),
      &cost_out, maxIter
    );
    if (res != 1)
      std::cerr << "WARNING: EMD_wrap did not converge\n";
    return cost_out;
}


double SolveOT_flat_full(
    const std::vector<double>& wx,
    const std::vector<double>& wy,
    const std::vector<std::vector<double>>& cost_matrix,
    const std::vector<int>& vx,
    const std::vector<int>& vy,
    int i0,                                            // base‑row offset
    int j0,                                            // base‑col offset
    const std::vector<std::vector<double>>& Vtplus,    // may be EMPTY_TABLE
    uint64_t maxIter = 100000
) {
    int n1 = wx.size(), n2 = wy.size(), sz = n1 * n2;
    D_buf.resize(sz); G_buf.resize(sz);
    α_buf.resize(n1); β_buf.resize(n2);

    // fill the flattened cost + future‐DP term
    for (int i = 0; i < n1; ++i) {
        int vi = vx[i];
        for (int j = 0; j < n2; ++j) {
            int vj = vy[j];
            D_buf[i*n2 + j] = cost_matrix[vi][vj]
                              + (!Vtplus.empty() ? Vtplus[i0 + i][j0 + j] : 0.0);
        }
    }

    // trivial 1×N or N×1
    if (n1 == 1 || n2 == 1) {
        double c = 0.0;
        for (int i = 0; i < n1; ++i)
            for (int j = 0; j < n2; ++j)
                c += wx[i] * D_buf[i*n2 + j] * wy[j];
        return c;
    }

    // otherwise EMD_wrap
    double cost_out = 0.0;
    int r = EMD_wrap(
      n1, n2,
      const_cast<double*>(wx.data()),
      const_cast<double*>(wy.data()),
      D_buf.data(), G_buf.data(),
      α_buf.data(), β_buf.data(),
      &cost_out, maxIter
    );
    if (r != 1) std::cerr<<"EMD_wrap failed\n";
    return cost_out;
}

// Empty table for the final step of the backward induction
static const std::vector<std::vector<double>> EMPTY_TABLE;


// REMOVED THE PRINTING HERE!!!
// Parameters:
//   X, Y         : Input paths (each row is a time step, columns are samples)
//   grid_size    : Grid size used for adapting/quantizing the paths.
//   markovian    : Switch between markovian (true) and full history (false) processing.
//   num_threads  : Number of threads to use (if <= 0, maximum available threads are used).
//   power        : Exponent for the cost function (only power 1 and 2 are optimized here).
double Nested(
    Eigen::MatrixXd& X,
    Eigen::MatrixXd& Y,
    double grid_size,
    const bool& markovian,
    int num_threads,
    const int power
) {
    // —————————————————————————————————————————————
    // 1) simulate/adapt/quantize exactly as before
    // —————————————————————————————————————————————
    int T        = X.rows() - 1;
    Eigen::MatrixXd adaptedX = path2adaptedpath(X, grid_size);
    Eigen::MatrixXd adaptedY = path2adaptedpath(Y, grid_size);

    // collect unique grid‑values
    std::set<double> v_set;
    v_set_add(adaptedX, v_set);
    v_set_add(adaptedY, v_set);

    // build quantization maps
    std::map<double,int> v2q;
    std::vector<double>   q2v;
    int pos = 0;
    for(double v: v_set){
        v2q[v] = pos++;
        q2v.push_back(v);
    }

    // quantize & lex‐sort
    Eigen::MatrixXi qX = sort_qpath(quantize_path(adaptedX, v2q).transpose());
    Eigen::MatrixXi qY = sort_qpath(quantize_path(adaptedY, v2q).transpose());

    // build the two conditional‐measure kernels
    auto mu_x = qpath2mu_x(qX, markovian);
    auto mu_y = qpath2mu_x(qY, markovian);
    auto kernel_x = mu_x2kernel_x(mu_x);
    auto kernel_y = mu_x2kernel_x(mu_y);

    // —————————————————————————————————————————————
    // 2) precompute base cost_matrix[i][j] = |q2v[i] - q2v[j]|^power
    // —————————————————————————————————————————————
    int V = (int)q2v.size();
    std::vector<std::vector<double>> cost_matrix(V, std::vector<double>(V));
    if(power == 1){
        for(int i=0;i<V;i++) for(int j=0;j<V;j++)
            cost_matrix[i][j] = std::abs(q2v[i] - q2v[j]);
    } else if(power == 2){
        for(int i=0;i<V;i++) for(int j=0;j<V;j++){
            double d = q2v[i] - q2v[j];
            cost_matrix[i][j] = d*d;
        }
    } else {
        for(int i=0;i<V;i++) for(int j=0;j<V;j++){
            cost_matrix[i][j] = std::pow(std::abs(q2v[i] - q2v[j]), (double)power);
        }
    }

    // —————————————————————————————————————————————
    // 3) allocate the value‐function table V[t][ix][iy]
    // —————————————————————————————————————————————
    std::vector<std::vector<std::vector<double>>> Vfunc(T);
    for(int t=0; t<T; t++){
        int nx = kernel_x[t].nc, ny = kernel_y[t].nc;
        Vfunc[t].assign(nx, std::vector<double>(ny, 0.0));
    }

    // —————————————————————————————————————————————
    // 4) set up OpenMP
    // —————————————————————————————————————————————
    if(num_threads <= 0) 
        num_threads = omp_get_max_threads();
    else
        num_threads = std::min(num_threads, omp_get_max_threads());
    omp_set_num_threads(num_threads);

    // —————————————————————————————————————————————
    // 5) backward‐induction in ONE parallel region
    // —————————————————————————————————————————————
    #pragma omp parallel
    {
      for (int t = T - 1; t >= 0; --t) {
        #pragma omp for collapse(2) schedule(dynamic)
        for (int ix = 0; ix < kernel_x[t].nc; ++ix) {
          for (int iy = 0; iy < kernel_y[t].nc; ++iy) {
            auto &dx = kernel_x[t].dists[ix];
            auto &dy = kernel_y[t].dists[iy];

            // If markovian: grab the next‐idx lists; else nullptr
            const std::vector<int>* x_idx = nullptr;
            const std::vector<int>* y_idx = nullptr;
            if (markovian && t < T - 1) {
              x_idx = &kernel_x[t].next_idx[ix];
              y_idx = &kernel_y[t].next_idx[iy];
            }

            // Choose which solver to call
            if (markovian) {
              // for final step, pass EMPTY_TABLE as Vnext
              const auto& Vnext = (t < T - 1 ? Vfunc[t+1] : EMPTY_TABLE);
              Vfunc[t][ix][iy] = SolveOT_flat(
                dx.weights, dy.weights,
                cost_matrix,
                dx.values, dy.values,
                x_idx, y_idx,
                Vnext
              );
            } else {
              // full‐history: compute base offsets i0,j0
              int i0 = kernel_x[t].nv_cums[ix];
              int j0 = kernel_y[t].nv_cums[iy];
              const auto& Vnext = (t < T - 1 ? Vfunc[t+1] : EMPTY_TABLE);
              Vfunc[t][ix][iy] = SolveOT_flat_full(
                dx.weights, dy.weights,
                cost_matrix,
                dx.values, dy.values,
                i0, j0,
                Vnext
              );
            }
          }
        }
      }
    } // omp parallel

    return Vfunc[0][0][0];
}
