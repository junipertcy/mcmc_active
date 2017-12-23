#include "Global.h"
#include <cmath>

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

template<typename VectorElemType, typename VectorIndexType>
static decltype(auto) qspartition(VectorElemType &A, VectorIndexType &B, int p, int r) {
    decltype(auto) x = A[r];
    decltype(auto) tempelem = A[0];  // TODO: better ways?
    decltype(auto) tempindex = B[0];
    int i = p - 1;
    for (int j = p; j <= r - 1; j++) {
        if (A[j] <= x) {
            i++;
            tempelem = A[i];
            tempindex = B[i];
            A[i] = A[j];
            B[i] = B[j];
            A[j] = tempelem;
            B[j] = tempindex;
        }
    }
    tempelem = A[i + 1];
    tempindex = B[i + 1];
    A[i + 1] = A[r];
    B[i + 1] = B[r];
    A[r] = tempelem;
    B[r] = tempindex;
    return i + 1;
};

template<typename VectorElemType, typename VectorIndexType>
static void quicksort(VectorElemType &A, VectorIndexType &B, int p, int r) {
    if (p < r) {
        int q = qspartition(A, B, p, r);
        quicksort(A, B, p, q - 1);
        quicksort(A, B, q + 1, r);
    }
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
