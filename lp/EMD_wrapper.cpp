/* This file is a c++ wrapper function for computing the transportation cost
 * between two vectors given a cost matrix.
 *
 * It was written by Antoine Rolet (2014) and mainly consists of a wrapper
 * of the code written by Nicolas Bonneel available on this page
 *          http://people.seas.harvard.edu/~nbonneel/FastTransport/
 *
 * It was then modified to make it more amenable to python inline calling
 *
 * Please give relevant credit to the original author (Nicolas Bonneel) if
 * you use this code for a publication.
 *
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
    std::cout << "Up to here" << std::endl;

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

    std::cout << ret << std::endl;
    std::cout << (int)net.OPTIMAL << std::endl;
    std::cout << (int)net.MAX_ITER_REACHED << std::endl;
    




    return ret;
}







int EMD_wrap_omp(int n1, int n2, double *X, double *Y, double *D, double *G,
             double* alpha, double* beta, double *cost, uint64_t maxIter, int numThreads)  {
    // beware M and C are stored in row major C style!!!

    using namespace lemon_omp;
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
    NetworkSimplexSimple<Digraph,double,double, node_id_type> net(di, true, (int) (n + m), n * m, maxIter, numThreads);

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


    return ret;
}


int main() {
    int n1 = 3, n2 = 3;
    double X[] = {0.3, 0.4, 0.3};  // Distribution 1 (supply)
    double Y[] = {0.2, 0.5, 0.3};  // Distribution 2 (demand)

    // Distance matrix D[n1][n2] in row-major order
    double D[] = {
        0.0, 1.0, 2.0,
        1.0, 0.0, 1.0,
        2.0, 1.0, 0.0
    };

    double G[9] = {0};            // Flow matrix (output)
    double alpha[3] = {0};          // Dual potentials (output)
    double beta[3] = {0};           // Dual potentials (output)
    double cost = 0.0;

    uint64_t maxIter = 100000;

    int result = EMD_wrap(n1, n2, X, Y, D, G, alpha, beta, &cost, maxIter);

    if (result == 1) {
        std::cout << "EMD cost: " << cost << std::endl;
        std::cout << "Flow matrix G:" << std::endl;
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                std::cout << G[i * n2 + j] << " ";
            }
            std::cout << std::endl;
        }
    } else {
        std::cerr << "EMD computation failed with codeeee: " << result << std::endl;
    }

    return 0;
}