#include "Global.h"
#include <cmath>
#include <random>

//[0,1)
static double rand01() {
    return ((double) rand() / ((double) RAND_MAX + 1));
};

//[0,n)
static int randN(int n) {
    return (int) ((double) rand() / ((double) RAND_MAX + 1) * n);
};

static double entropy(double_vec_t dist, unsigned int size) {
    double temp = 0.0;
    for (unsigned int i = 0; i < size; ++i) {
        if (dist[i] > 0) {
            temp -= dist[i] * log(dist[i]);// log: natural logarithm
        }
    }
    return temp;
};

static int getIndexProb(double_vec_t &probarray, int dimension) {
    int m;
    int l = 0;
    int r = dimension - 1;
    double dRandValue = (rand01()) * (probarray[r]);
    do {
        m = (l + r) / 2;
        if (dRandValue < probarray[m]) r = m;
        if (dRandValue >= probarray[m]) l = m + 1;
    } while (l != r);
    return l;
};
