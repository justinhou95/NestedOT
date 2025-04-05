#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include <Eigen/Dense>
#include <set>
#include <map>
#include <algorithm>
#include "header.h"

Eigen::MatrixXd Lmatrix2paths(
    Eigen::MatrixXd L,
    int n_sample,
    bool normalize = false,
    int seed = 0,
    bool verbose = false
) {

    int T = L.rows();

    // Random noise matrix (T, n_sample)
    std::mt19937 gen(seed);
    std::normal_distribution<> d(0.0, 1.0);
    Eigen::MatrixXd noise1(T, n_sample);
    for (int i = 0; i < T; ++i)
        for (int j = 0; j < n_sample; ++j)
            noise1(i, j) = d(gen);

    // X = L * noise1
    Eigen::MatrixXd X = L * noise1;

    // Add a row of zeros at the top (T+1, n_sample)
    Eigen::MatrixXd X_extended(T + 1, n_sample);
    X_extended.row(0).setZero();
    X_extended.block(1, 0, T, n_sample) = X;

    return X_extended;
}

Eigen::MatrixXd path2adaptedpath(const Eigen::MatrixXd& X, double delta_n) {
    Eigen::MatrixXd adaptedX = X;
    for (int i = 0; i < adaptedX.rows(); ++i) {
        for (int j = 0; j < adaptedX.cols(); ++j) {
            adaptedX(i, j) = std::floor(X(i, j) / delta_n + 0.5) * delta_n;
        }
    }
    return adaptedX;
}


void v_set_add(const Eigen::MatrixXd& mat, std::set<double>& unique_set) {
    for (int i = 0; i < mat.rows(); ++i) {
        for (int j = 0; j < mat.cols(); ++j) {
            unique_set.insert(mat(i, j));
        }
    }
}

Eigen::MatrixXi quantize_path(Eigen::MatrixXd& adaptedX, std::map<double, int>& v2q) {
    Eigen::MatrixXi quantizedX(adaptedX.rows(), adaptedX.cols());

    for (int i = 0; i < adaptedX.rows(); ++i) {
        for (int j = 0; j < adaptedX.cols(); ++j) {
            double v = adaptedX(i,j);
            quantizedX(i,j) = v2q[v];
        }
    }
    return quantizedX;
}


Eigen::MatrixXi sort_qpath(const Eigen::MatrixXi& path) {
    int n_rows = path.rows();
    int n_cols = path.cols();

    std::vector<int> indices(n_rows);
    std::iota(indices.begin(), indices.end(), 0);  // Fill with 0 to n_rows-1

    auto row_cmp = [&](int i, int j) {
        for (int col = 0; col < n_cols; ++col) {
            if (path(i, col) < path(j, col)) return true;
            if (path(i, col) > path(j, col)) return false;
        }
        return false; // rows are equal
    };

    std::sort(indices.begin(), indices.end(), row_cmp);

    Eigen::MatrixXi sorted(n_rows, n_cols);
    for (int i = 0; i < n_rows; ++i) {
        sorted.row(i) = path.row(indices[i]);
    }

    return sorted;
}

std::vector<std::map<std::vector<int>, std::map<int, int>>> qpath2mu_x(Eigen::MatrixXi& qpath) {
    int T = qpath.cols()-1;

    std::vector<std::map<std::vector<int>, std::map<int, int>>> mu_x(T);

    for (int t = 0; t < T; ++t) {
        for (int i = 0; i < qpath.rows(); ++i) {
            std::vector<int> pre_path;
            for (int k = 0; k <= t; ++k) {
                pre_path.push_back(qpath(i, k));
            }
            int next_val = qpath(i, t + 1);
            mu_x[t][pre_path][next_val] += 1;
        }
    }
    return mu_x;
}

std::vector<ConditionalDistribution> mu_x2kernel_x(std::vector<std::map<std::vector<int>, std::map<int, int>>>& mu_x){
    int T =  mu_x.size();
    std::vector<ConditionalDistribution> kernel_x(T);
    for (int t = 0; t < T; t++){
        std::map<std::vector<int>, std::map<int, int>>& mu_x_t = mu_x[t];
        kernel_x[t].nc = mu_x_t.size();
        kernel_x[t].t = t;

        kernel_x[t].nv_cums.push_back(0);
        for (auto pair : mu_x_t) {
            const std::vector<int>& condition = pair.first;
            const std::map<int, int>& distribution = pair.second;

            std::vector<int> values;
            std::vector<int> counts;
            std::vector<double> weights(distribution.size());
            int sum = 0;
            for (auto d : distribution){
                values.push_back(d.first);
                counts.push_back(d.second);
                sum += d.second;
            }
            for (int i=0; i < weights.size(); i++){
                weights[i] = (double) counts[i]/sum;                
            }

            const Distribution dist = {values, weights}; 
            kernel_x[t].conds.push_back(condition);
            kernel_x[t].dists.push_back(dist);
            kernel_x[t].nvs.push_back(distribution.size());
            kernel_x[t].nv_cums.push_back(kernel_x[t].nv_cums.back() + distribution.size());
        }
    }
    return kernel_x;
}