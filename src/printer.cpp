#include <iostream>
#include <vector>
#include <map>
#include "header_dist.h"


// Helper function to print a vector<int>
void print_vector_int(const std::vector<int>& vec) {
    std::cout << "(";
    for (int i = 0; i < vec.size(); ++i) {
        std::cout << vec[i] << ",";
    }
    std::cout << ")";
}

void print_vector_double(const std::vector<double>& vec) {
    std::cout << "(";
    for (int i = 0; i < vec.size(); ++i) {
        std::cout << vec[i] << ",";
    }
    std::cout << ")";
}

void print_map(const std::map<int, int>& map) {
    std::cout << "{";
    for (const auto& inner_pair : map) {
        std::cout << inner_pair.first << ":" << inner_pair.second << ",";
    }
    std::cout << "}\n";
}

void print_mu_x(const std::map<std::vector<int>, std::map<int, int>>& mu_x) {
    for (const auto& outer_pair : mu_x) {
        print_vector_int(outer_pair.first);  // print the key vector
        std::cout << " -> ";
        print_map(outer_pair.second);
    }
}

// Print condition as tuple
void printCondition(std::vector<int> cond){
    std::cout << "Condition: ";
    print_vector_int(cond);
    std::cout << " ";
}

// Print Values and Weight of Distribution as lists
void printDistribution(Distribution dist){
    std::cout << "Values: ";
    print_vector_int(dist.values);
    std::cout << " Weights:";
    print_vector_double(dist.weights);
    std::cout << std::endl;
}


void print_kernel_x(std::vector<ConditionalDistribution>& kernel_x){
    int T = kernel_x.size();
    for(int t = 0; t < T; t++){
        std::cout << t << std::endl;
        std::cout << kernel_x[t].nc << std::endl;
        // print_vector_int(kernel_x[t].nvs);
        // std::cout << std::endl;
        // print_vector_int(kernel_x[t].nv_cums);
        std::cout << std::endl;
        print_map(kernel_x[t].conds2idx);
        std::cout << std::endl;
        if (t < T-1){
            for (std::vector<int> idx_list : kernel_x[t].next_idx){
                print_vector_int(idx_list);
                std::cout << std::endl;
            }
        }
        std::cout << std::endl;

        for (int ix =0; ix < kernel_x[t].nc; ix++){
            printCondition(kernel_x[t].conds[ix]);
            printDistribution(kernel_x[t].dists[ix]);
            std::cout << std::endl;
        }


    }
}