#pragma once

#include <vector>

struct Distribution{
    std::vector<int> values;
    std::vector<float> weights;
};

struct ConditionalDistribution{
    int t;
    int nc;
    std::vector<std::vector<int>> conds;
    std::vector<Distribution> dists;
    std::vector<int> nvs;
    std::vector<int> nv_cums;
};