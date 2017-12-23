#include "Learner.h"
#include "Utility.cpp"

const unsigned int Learner::QUERYSTRATEGY = 1;

Learner::Learner(const MCMC &mca, const MCMC &mcb, const set<unsigned int> &topVtxSet):
        m_MC_A(mca), m_MC_B(mcb), m_topVtxSet(topVtxSet)
{
    m_numVtx = mca.getTypeModel().getGraph().getNumVtx();
    m_numType = mca.getTypeModel().getNumType();

    m_arrayTopVtxNo.resize(m_numVtx);
    m_arrayVtxNo.resize(m_numVtx);
    m_arrayLearnScores.resize(m_numVtx);

    m_arrayVtxNoSort.resize(m_numVtx);
    m_arrayLearnScoresSort.resize(m_numVtx);
}

unsigned int MutualInfo::getNumVtx() const noexcept { return Learner::getNumVtx(); }

unsigned int MutualInfo::getNumTypeGraph() const noexcept { return Learner::getNumTypeGraph(); }

unsigned int MutualInfo::getNumTypeModel() const noexcept { return Learner::getNumTypeModel(); }

void MutualInfo::resetForNextInit() noexcept { ; }

MutualInfo::MutualInfo(const MCMC &mca, const MCMC &mcb, const set<unsigned> &topVtxSet) : Learner(mca, mcb, topVtxSet) {
    m_accumuCondEntropy.resize(m_numVtx);
    m_numAccumuCondEntropy.resize(m_numVtx);
    m_accumuMargDistri.resize(m_numVtx);

    for (unsigned int i = 0; i < m_numVtx; ++i) {
        m_accumuMargDistri[i].resize(m_numType);  // number of groups
    }
}


void MutualInfo::resetForNextPhase() noexcept {

    unsigned i, j;
    for (i = 0; i < m_numVtx; i++) {
        for (j = 0; j < m_numType; j++) {
            m_accumuMargDistri[i][j] = 0.0;
        }
    }
    for (i = 0; i < m_numVtx; i++) {
        m_accumuCondEntropy[i] = 0.0;
        m_numAccumuCondEntropy[i] = 0;
    }
}

unsigned MutualInfo::getTopVtx(uint_vec_t &arrayTop, unsigned numTop) noexcept {
    unsigned i, j;
    double margEntropy, condEntropy, mutualEntropy;
    for (i = 0, m_arraysize = 0; i < m_numVtx; i++) {
        if (m_topVtxSet.count(i)) {
            continue;
        }
        if (m_numAccumuCondEntropy[i] == 0) {
            cerr << "need more MCMC steps" << endl;
            continue;
        }
        condEntropy = m_accumuCondEntropy[i] / (double) m_numAccumuCondEntropy[i];

        // To normalize the accumulated Marginal Distribution
        double dSum = 0.0;
        for (j = 0; j < m_numType; j++) {
            dSum += m_accumuMargDistri[i][j];
        }
        if (dSum == 0.0) {
            cerr << "need more MCMC steps" << endl;
            continue;
        }
        for (j = 0; j < m_numType; j++) {
            m_accumuMargDistri[i][j] /= dSum;
        }
        margEntropy = entropy(m_accumuMargDistri[i], m_numType);  // Now, we can directly compute the average entropy

        mutualEntropy = margEntropy - condEntropy;
        m_arrayLearnScores[m_arraysize] = mutualEntropy;
        m_arrayLearnScoresSort[m_arraysize] = mutualEntropy;
        m_arrayVtxNo[m_arraysize] = i;
        m_arrayVtxNoSort[m_arraysize++] = i;
    }
    quicksort<double_vec_t, uint_vec_t>(m_arrayLearnScoresSort, m_arrayVtxNoSort, 0, m_arraysize - 1);
    m_numtop = 0;
    for (i = m_arraysize - 1; m_numtop < numTop; --i) {
        m_arrayTopVtxNo[m_numtop++] = m_arrayVtxNoSort[i];
    }
    cout << "tops:" << endl;
    for (i = 0; i < m_numtop; i++) {
        arrayTop[i] = m_arrayTopVtxNo[i];
        cout << arrayTop[i] << '\t';
    }
    cout << endl;
    return m_numtop;
}

void MutualInfo::updateData() noexcept {
    unsigned mutatevtxno_a = m_MC_A.getMutateVtx();
    unsigned mutatevtxno_b = m_MC_B.getMutateVtx();
    const vector <pair<unsigned int, double>> &lhvpairs_a = m_MC_A.getLHVariPairs();
    const vector <pair<unsigned int, double>> &lhvpairs_b = m_MC_B.getLHVariPairs();

    double_vec_t probarray;
    probarray.resize(lhvpairs_a.size());

    for (unsigned int i = 0; i < lhvpairs_a.size(); i++) {
        m_accumuMargDistri[mutatevtxno_a][lhvpairs_a[i].first] += lhvpairs_a[i].second;
    }
    for (unsigned int i = 0; i < lhvpairs_b.size(); i++) {
        m_accumuMargDistri[mutatevtxno_b][lhvpairs_b[i].first] += lhvpairs_b[i].second;
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
