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

void print_mu_x_markovian(const std::map<int, std::map<int, int>>& mu_x) {
    for (const auto& outer_pair : mu_x) {
        std::cout << outer_pair.first << " -> ";
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

