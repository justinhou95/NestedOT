#pragma once

#include <iostream>
#include <vector>
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
    std::vector<int> nv_cums; // Only use this for non-markovian
    
    std::map<int,int> v2idx; // Only use this for markovian
    std::vector<std::vector<int>> next_idx; 
};