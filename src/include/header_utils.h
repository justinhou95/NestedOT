#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <set>
#include <map>
#include <random>
#include <algorithm>

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