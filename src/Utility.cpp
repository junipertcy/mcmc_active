#include "Global.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <set>

static const double LARGE_NEGATIVE_VALUE = -(numeric_limits<double>::max());

//[0,1)
static double rand01() {
    return ((double) rand() / ((double) RAND_MAX + 1));
};

//[0,n)
static int randN(int n) {
    return (int) ((double) rand() / ((double) RAND_MAX + 1) * n);
};

static ostream &print_ui_set(const set<unsigned> &uiset, ostream &os = cout) {
    for (auto iter = uiset.begin(); iter != uiset.end(); iter++) {
        os << *iter << endl;
    }
    return os;
};

static ostream &print_ui_vec(const vector<unsigned> &uivec, ostream &os = cout) {
    unsigned int sum = 0;
    for (auto iter = uivec.begin(); iter != uivec.end(); iter++) {
        os << *iter << " ";
        sum += *iter;
    }
    os << endl;
    os << "Average value of the vector is: " << (double) sum / uivec.size() << endl;
    return os;
};


static double entropy(double_vec_t dist, unsigned size) {
    double temp = 0.0;
    for (unsigned int i = 0; i < size; ++i) {
        if (dist[i] > 0) {
            temp -= dist[i] * log(dist[i]);// log: natural logarithm
        }
    }
    return temp;
};

static void getTopN(const double *array, int *top, int size, int n) {
    int i;
    int topnum = n;
    if (size < n) {
        cerr << "finding top n: n is over the array size boundary" << endl;
        topnum = size;
    }
    auto *localArray = new double[size];
    for (i = 0; i < size; i++)
        localArray[i] = array[i];
    int temp;
    for (i = 0; i < topnum; i++) {
        temp = 0;
        for (int j = 1; j < size; j++) {
            if (localArray[j] > localArray[temp])
                temp = j;
        }
        top[i] = temp;
        localArray[temp] = LARGE_NEGATIVE_VALUE; //mark as already picked
    }
    delete[] localArray;
};

template<typename ArrayElemType, typename ArrayIndexType>
static int qspartition(ArrayElemType *A, ArrayIndexType *B, int p, int r) {
    ArrayElemType x = A[r];
    ArrayElemType tempelem;
    ArrayIndexType tempindex;
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

template<typename ArrayElemType, typename ArrayIndexType>
static void quicksort(ArrayElemType *A, ArrayIndexType *B, int p, int r) {
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

template<typename fstElemType, typename sedElemType>
class Larger_Pair_Second {
public:
    bool operator()(const pair<fstElemType, sedElemType> &p1, const pair<fstElemType, sedElemType> &p2) {
        return p1.second > p2.second;
    }
};

template<class T>
static std::string to_string(const T &t) {
    std::stringstream ss;
    ss << t;
    return ss.str();
};

static unsigned sizeCommon(vector<unsigned> &vec1, vector<unsigned> &vec2) {
    unsigned i, j;
    unsigned sizecommon = 0;
    for (i = 0; i < vec1.size(); i++) {
        for (j = 0; j < vec2.size(); j++) {
            if (vec1[i] == vec2[j]) {
                sizecommon++;
                break;
            }
        }
    }
    return sizecommon;
};

static double
calcCondEntropy(unsigned NUM_TYPE_C1, unsigned NUM_TYPE_C2, unsigned NUM_VTX, unsigned *C1, unsigned *C2) {
    unsigned int i, j;
    double H_C1, H_C2, MI, VI, CVI, HC1CondC2, HC1CondC2Norm;
    //
    vector<vector<unsigned int>> clustering1;
    vector<vector<unsigned int>> clustering2;
    //
    for (i = 0; i < NUM_TYPE_C1; i++) {
        vector<unsigned int> cluster;
        clustering1.push_back(cluster);
    }
    for (i = 0; i < NUM_TYPE_C2; i++) {
        vector<unsigned int> cluster;
        clustering2.push_back(cluster);
    }
    //
    for (i = 0; i < NUM_VTX; i++) {
        clustering1[C1[i]].push_back(i);
        clustering2[C2[i]].push_back(i);
    }
    //
    double_vec_t distri_c1;
    distri_c1.resize(NUM_TYPE_C1);

    double_vec_t distri_c2;
    distri_c2.resize(NUM_TYPE_C2);

    for (i = 0; i < NUM_TYPE_C1; i++) {
        distri_c1[i] = float(clustering1[i].size()) / NUM_VTX;
    }
    H_C1 = entropy(distri_c1, NUM_TYPE_C1);
    for (i = 0; i < NUM_TYPE_C2; i++) {
        distri_c2[i] = float(clustering2[i].size()) / NUM_VTX;
    }
    H_C2 = entropy(distri_c2, NUM_TYPE_C2);
    MI = 0.0;
    double pij;
    double pi;
    double pj;
    for (i = 0; i < NUM_TYPE_C1; i++) {
        for (j = 0; j < NUM_TYPE_C2; j++) {
            pij = double(sizeCommon(clustering1[i], clustering2[j])) / NUM_VTX;
            if (pij == 0)
                continue;
            pi = double(clustering1[i].size()) / NUM_VTX;
            pj = double(clustering2[j].size()) / NUM_VTX;
            MI += pij * log(pij / pi / pj);
        }
    }
    VI = H_C1 + H_C2 - 2 * MI;
    HC1CondC2 = H_C1 - MI;
    HC1CondC2Norm = (H_C1 - MI) / H_C1;
    CVI = VI / log(double(NUM_VTX));
    return HC1CondC2Norm;
};

static double calcNormMI(unsigned NUM_TYPE_C1, unsigned NUM_TYPE_C2, unsigned NUM_VTX, unsigned *C1, unsigned *C2) {
    unsigned i, j;
    double H_C1, H_C2, MI;
    //
    vector<vector<unsigned int>> clustering1;
    vector<vector<unsigned int>> clustering2;
    //
    for (i = 0; i < NUM_TYPE_C1; i++) {
        vector<unsigned int> cluster;
        clustering1.push_back(cluster);
    }
    for (i = 0; i < NUM_TYPE_C2; i++) {
        vector<unsigned int> cluster;
        clustering2.push_back(cluster);
    }
    //
    for (i = 0; i < NUM_VTX; i++) {
        clustering1[C1[i]].push_back(i);
        clustering2[C2[i]].push_back(i);
    }
    //
    double_vec_t distri_c1;
    double_vec_t distri_c2;

    distri_c1.resize(NUM_TYPE_C1);
    distri_c2.resize(NUM_TYPE_C2);

    for (i = 0; i < NUM_TYPE_C1; i++) {
        distri_c1[i] = float(clustering1[i].size()) / NUM_VTX;
    }
    H_C1 = entropy(distri_c1, NUM_TYPE_C1);
    for (i = 0; i < NUM_TYPE_C2; i++) {
        distri_c2[i] = float(clustering2[i].size()) / NUM_VTX;
    }
    H_C2 = entropy(distri_c2, NUM_TYPE_C2);
    MI = 0.0;
    double pij;
    double pi;
    double pj;
    for (i = 0; i < NUM_TYPE_C1; i++) {
        for (j = 0; j < NUM_TYPE_C2; j++) {
            pij = double(sizeCommon(clustering1[i], clustering2[j])) / NUM_VTX;
            if (pij == 0)
                continue;
            pi = double(clustering1[i].size()) / NUM_VTX;
            pj = double(clustering2[j].size()) / NUM_VTX;
            MI += pij * log(pij / pi / pj);
        }
    }
    return MI / H_C1;
};
