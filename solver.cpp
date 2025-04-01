#include <iostream>
#include <vector>
#include "header.h"

std::vector<std::vector<float>> SquareCost(std::vector<int>& vx, std::vector<int>& vy, std::vector<std::vector<float>>& cost_matrix){
    int m = vx.size();
    int n = vy.size();
    std::vector<std::vector<float>> cost(m, std::vector<float>(n, 0.0f));
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            cost[i][j] = cost_matrix[vx[i]][vy[j]];
        }
    }
    return cost;
}

void AddDppValue(std::vector<std::vector<float>>& cost, std::vector<std::vector<float>>& Vtplus, int& i0, int& j0){
    int m = cost.size();
    int n = cost[0].size();
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            cost[i][j] += Vtplus[i0+i][j0+j];
        }
    }
}

float SolveOT(std::vector<float>& wx, std::vector<float>& wy, std::vector<std::vector<float>>& cost){
    int n1 = wx.size(); 
    int n2 = wy.size();

    float* X = wx.data();  
    float* Y = wy.data();

    std::vector<float> flatcost(cost.size() * cost[0].size());
    for (const auto& row : cost) {
        flatcost.insert(flatcost.end(), row.begin(), row.end());
    }
    float* D = flatcost.data();

    std::vector<float> Gvec(n1*n2);
    float* G = Gvec.data();            // Flow matrix (output)
    std::vector<float> alphavec(n1);
    float* alpha = alphavec.data();          // Dual potentials (output)
    std::vector<float> betavec(n2);
    float* beta = betavec.data();           // Dual potentials (output)
    float c = 0.0;

    uint64_t maxIter = 100000;

    // int result = EMD_wrap(n1, n2, X, Y, D, G, alpha, beta, &c, maxIter);
    return c;
}


