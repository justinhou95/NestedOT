#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <set>
#include <map>

struct Distribution{
    std::vector<int> values;
    std::vector<double> weights;
};

struct ConditionalDistribution{
    int t;
    int nc;
    std::vector<std::vector<int>> conds;
    std::vector<Distribution> dists;
    std::vector<int> nvs;
    std::vector<int> nv_cums;
};


Eigen::MatrixXd path2adaptedpath(const Eigen::MatrixXd& samples, double delta_n);

void v_set_add(const Eigen::MatrixXd& mat, std::set<double>& unique_set);

Eigen::MatrixXi quantize_path(Eigen::MatrixXd& adaptedX, std::map<double, int>& v2q);

Eigen::MatrixXi sort_qpath(const Eigen::MatrixXi& path);

std::vector<std::map<std::vector<int>, std::map<int, int>>> qpath2mu_x(Eigen::MatrixXi& qpath);

std::vector<ConditionalDistribution> mu_x2kernel_x(std::vector<std::map<std::vector<int>, std::map<int, int>>>& mu_x);