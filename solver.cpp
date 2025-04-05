#include <iostream>
#include <vector>
#include "header.h"

/* 
 * Please give relevant credit to the original author (Nicolas Bonneel) if
 * you use this code for a publication.
 */


 #include "network_simplex_simple.h"
 #include "network_simplex_simple_omp.h"
 #include "EMD.h"
 #include <cstdint>
 
 
 int EMD_wrap(int n1, int n2, double *X, double *Y, double *D, double *G,
                 double* alpha, double* beta, double *cost, uint64_t maxIter)  {
     // beware M and C are stored in row major C style!!!
 
     using namespace lemon;
     uint64_t n, m, cur;
 
     typedef FullBipartiteDigraph Digraph;
     DIGRAPH_TYPEDEFS(Digraph);
 
     // Get the number of non zero coordinates for r and c
     n=0;
     for (int i=0; i<n1; i++) {
         double val=*(X+i);
         if (val>0) {
             n++;
         }else if(val<0){
             return INFEASIBLE;
         }
     }
     m=0;
     for (int i=0; i<n2; i++) {
         double val=*(Y+i);
         if (val>0) {
             m++;
         }else if(val<0){
             return INFEASIBLE;
         }
     }
 
     // Define the graph
 
     std::vector<uint64_t> indI(n), indJ(m);
     std::vector<double> weights1(n), weights2(m);
     Digraph di(n, m);
     NetworkSimplexSimple<Digraph,double,double, node_id_type> net(di, true, (int) (n + m), n * m, maxIter);
 
     // Set supply and demand, don't account for 0 values (faster)
 
     cur=0;
     for (uint64_t i=0; i<n1; i++) {
         double val=*(X+i);
         if (val>0) {
             weights1[ cur ] = val;
             indI[cur++]=i;
         }
     }
 
     // Demand is actually negative supply...
 
     cur=0;
     for (uint64_t i=0; i<n2; i++) {
         double val=*(Y+i);
         if (val>0) {
             weights2[ cur ] = -val;
             indJ[cur++]=i;
         }
     }
 
 
     net.supplyMap(&weights1[0], (int) n, &weights2[0], (int) m);
 
     // Set the cost of each edge
     int64_t idarc = 0;
     for (uint64_t i=0; i<n; i++) {
         for (uint64_t j=0; j<m; j++) {
             double val=*(D+indI[i]*n2+indJ[j]);
             net.setCost(di.arcFromId(idarc), val);
             ++idarc;
         }
     }
 
     // Solve the problem with the network simplex algorithm
 
     int ret=net.run();
     uint64_t i, j;
     if (ret==(int)net.OPTIMAL || ret==(int)net.MAX_ITER_REACHED) {
         *cost = 0;
         Arc a; di.first(a);
         for (; a != INVALID; di.next(a)) {
             i = di.source(a);
             j = di.target(a);
             double flow = net.flow(a);
             *cost += flow * (*(D+indI[i]*n2+indJ[j-n]));
             *(G+indI[i]*n2+indJ[j-n]) = flow;
             *(alpha + indI[i]) = -net.potential(i);
             *(beta + indJ[j-n]) = net.potential(j);
         }
 
     }
     
    //  std::cout << "Some information" << std::endl;
    //  std::cout << ret << std::endl;
    //  std::cout << (int)net.OPTIMAL << std::endl;
    //  std::cout << (int)net.MAX_ITER_REACHED << std::endl;
     
     return ret;
 }

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

double SolveOT(std::vector<double>& wx, std::vector<double>& wy, std::vector<std::vector<double>>& cost){
    // We need to cast the double to double or have all setup in double in the very beginning, I would prefer the later option for accuracy.
    int n1 = wx.size(); 
    int n2 = wy.size();

    double* X = wx.data();  
    double* Y = wy.data();

    std::vector<double> flatcost(cost.size() * cost[0].size());
    for (const auto& row : cost) {
        flatcost.insert(flatcost.end(), row.begin(), row.end());  // Is it necessary to do like this or we can directly pass the vector double data to a double pointer
    }
    double* D = flatcost.data();

    std::vector<double> Gvec(n1*n2);
    double* G = Gvec.data();            // Flow matrix (output)
    std::vector<double> alphavec(n1);
    double* alpha = alphavec.data();          // Dual potentials (output)
    std::vector<double> betavec(n2);
    double* beta = betavec.data();           // Dual potentials (output)
    double c = 0.0;

    uint64_t maxIter = 100000;

    int result = EMD_wrap(n1, n2, X, Y, D, G, alpha, beta, &c, maxIter);
    return c;
}

 

