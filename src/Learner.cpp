#include "Learner.h"
#include "output_functions.h"
#include "Utility.cpp"


Learner::Learner(const MCMC &mca, const MCMC &mcb, const set<unsigned int> &topVtxSet):
        m_MC_A(mca),
        m_MC_B(mcb),
        m_topVtxSet(topVtxSet)
{
//    m_MC_A = std::make_unique<MCMC>(*mca);
//    m_MC_B = *mcb;

    N_ = mca.getTypeModel().getGraph().getNumVtx();
    Q_ = mca.getTypeModel().getNumType();


}

unsigned int MutualInfo::getNumVtx() const noexcept { return Learner::getNumVtx(); }

unsigned int MutualInfo::getNumTypeGraph() const noexcept { return Learner::getNumTypeGraph(); }

unsigned int MutualInfo::getNumTypeModel() const noexcept { return Learner::getNumTypeModel(); }

void MutualInfo::resetForNextInit() noexcept { ; }

MutualInfo::MutualInfo(const MCMC &mca, const MCMC &mcb, const set<unsigned> &topVtxSet) : Learner(mca, mcb, topVtxSet) {
    m_accumuCondEntropy.resize(N_);
    m_numAccumuCondEntropy.resize(N_);
    m_accumuMargDistri.resize(N_);

    for (unsigned int i = 0; i < N_; ++i) {
        m_accumuMargDistri[i].resize(Q_);  // number of groups
    }
}


void MutualInfo::resetForNextPhase() noexcept {
    for (unsigned int i = 0; i < N_; i++) {
        for (unsigned int j = 0; j < Q_; j++) {
            m_accumuMargDistri[i][j] = 0.0;
        }
    }
    for (unsigned int i = 0; i < N_; i++) {
        m_accumuCondEntropy[i] = 0.0;
        m_numAccumuCondEntropy[i] = 0;
    }
}

unsigned MutualInfo::getTopVtx(uint_vec_t &arrayTop, unsigned int numTop) noexcept {
    unsigned i, j, m_arraysize;
    double margEntropy, condEntropy, mutualEntropy;

    uint_vec_t m_arrayTopVtxNo; //some query strategies may output more than "numTop" vertices
    m_arrayTopVtxNo.resize(N_);

    uint_vec_t m_arrayVtxNoSort;
    m_arrayVtxNoSort.resize(N_);

    double_vec_t m_arrayLearnScoresSort;
    m_arrayLearnScoresSort.resize(N_);

    for (i = 0, m_arraysize = 0; i < N_; ++i) {
        if (m_topVtxSet.count(i)) {
            continue;
        }
        if (m_numAccumuCondEntropy[i] == 0) {
            cerr << "need more MCMC steps" << endl;
            continue;
        }
        condEntropy = m_accumuCondEntropy[i] / (double) m_numAccumuCondEntropy[i];
//        output_vec<float_vec_t>(m_numAccumuCondEntropy);
        // To normalize the accumulated Marginal Distribution
        double dSum = 0.0;
        for (j = 0; j < Q_; ++j) {
            dSum += m_accumuMargDistri[i][j];
        }
        if (dSum == 0.0) {
            cerr << "need more MCMC steps" << endl;
            continue;
        }
        for (j = 0; j < Q_; ++j) {
            m_accumuMargDistri[i][j] /= dSum;
        }
        margEntropy = entropy(m_accumuMargDistri[i], Q_);  // Now, we can directly compute the average entropy

        mutualEntropy = margEntropy - condEntropy;
        m_arrayLearnScoresSort[m_arraysize] = mutualEntropy;
        m_arrayVtxNoSort[m_arraysize++] = i;
    }

    quicksort<double_vec_t, uint_vec_t>(m_arrayLearnScoresSort, m_arrayVtxNoSort, 0, m_arraysize - 1);
//    std::clog << "----\n";
//    std::clog << "The gain in mutual information (the more the better) of each node is: \n";
    unsigned int m_numtop = 0;
    for (i = m_arraysize - 1; m_numtop < numTop; --i) {
        m_arrayTopVtxNo[m_numtop++] = m_arrayVtxNoSort[i];
    }

    std::clog << "tops:" << endl;
    for (i = 0; i < m_numtop; i++) {
        arrayTop[i] = m_arrayTopVtxNo[i];
        clog << arrayTop[i] << '\t';
    }
    std::clog << endl;
    return m_numtop;
}

void MutualInfo::updateData() noexcept {
    unsigned int mutatevtxno_a = m_MC_A.getMutateVtx();
    unsigned int mutatevtxno_b = m_MC_B.getMutateVtx();
    const vector <pair<unsigned int, double>> &lhvpairs_a = m_MC_A.getLHVariPairs();
    const vector <pair<unsigned int, double>> &lhvpairs_b = m_MC_B.getLHVariPairs();

    double_vec_t probarray;
    probarray.resize(lhvpairs_a.size());

    for (auto la: lhvpairs_a) {
        m_accumuMargDistri[mutatevtxno_a][la.first] += la.second;
    }

    for (auto lb: lhvpairs_b) {
        m_accumuMargDistri[mutatevtxno_b][lb.first] += lb.second;
    }

    for (unsigned int i = 0; i < lhvpairs_a.size(); i++) {
        probarray[i] = lhvpairs_a[i].second;
    }

    m_accumuCondEntropy[mutatevtxno_a] += entropy(probarray, (unsigned) (lhvpairs_a.size()));
    m_numAccumuCondEntropy[mutatevtxno_a]++;

    for (unsigned int i = 0; i < lhvpairs_b.size(); i++) {
        probarray[i] = lhvpairs_b[i].second;
    }
    m_accumuCondEntropy[mutatevtxno_b] += entropy(probarray, (unsigned) (lhvpairs_b.size()));
    m_numAccumuCondEntropy[mutatevtxno_b]++;

}
