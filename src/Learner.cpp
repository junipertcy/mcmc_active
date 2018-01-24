#include <algorithm> // For sort()
#include "Learner.h"
#include "output_functions.h"
#include "Utility.cpp"

bool compare(const pair<double,unsigned int>&i, const pair<double,unsigned int>&j) noexcept {
    return i.first > j.first;
}

Learner::Learner(const MCMC &mca, const MCMC &mcb, const set<unsigned int> &topVtxSet):
        m_MC_A(mca),
        m_MC_B(mcb),  // TODO: why!?
        m_topVtxSet(topVtxSet) {
    N_ = mca.getTypeModel().getGraph().get_N();
    Q_ = mca.getTypeModel().get_Q();
}

unsigned int MutualInfo::get_N() const noexcept { return Learner::get_N(); }

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
    unsigned i, j;
    double margEntropy, condEntropy, mutualEntropyGain;
    vector <pair<double,unsigned int> > nodes_mi_gains;

    nodes_mi_gains.resize(N_);
    for (i = 0; i < N_; ++i) {
        if (m_topVtxSet.count(i)) {
            nodes_mi_gains[i] = pair<double,unsigned int>(0., i);
            continue;
        }
        if (m_numAccumuCondEntropy[i] == 0) {
            cerr << "need more MCMC steps" << endl;
            continue;
        }
        condEntropy = m_accumuCondEntropy[i] / (double) m_numAccumuCondEntropy[i];

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
        // MI(v) = H(v) − H(v | G \ v)
        //       = entropy of the average conditional distribution (margEntropy) - average of the conditional entropy (condEntropy);

        mutualEntropyGain = margEntropy - condEntropy;
//        std::clog << "[DEBUG] @i=" << i <<"; margEntropy: " << margEntropy << "; condEntropy: " << condEntropy << "\n";

        nodes_mi_gains[i] = pair<double,unsigned int>(mutualEntropyGain, i);
    }

    for (unsigned int node_id = 0; node_id < N_; ++node_id) {
        if (node_id > 0) {
            std::cout << ",";
        }
        std::cout << nodes_mi_gains[node_id].first;
    }
    std::cout << "\n";

    std::sort(nodes_mi_gains.begin(),nodes_mi_gains.end(), compare);

    std::clog << "The indexes of the selected node(s) are: ";
    for (i = 0; i < numTop; ++i) {
        arrayTop[i] = nodes_mi_gains.at(i).second;
        clog << arrayTop[i] << ',';
    }
    std::clog << "\n";

    return numTop;
}

void MutualInfo::updateData() noexcept {
    unsigned int selected_vtx_a = m_MC_A.get_selected_vtx();
    unsigned int selected_vtx_b = m_MC_B.get_selected_vtx();
    const vector <pair<unsigned int, double>> &lhvpairs_a = m_MC_A.get_likelihood_variation_pairs();
    const vector <pair<unsigned int, double>> &lhvpairs_b = m_MC_B.get_likelihood_variation_pairs();

    double_vec_t probarray;
    probarray.resize(lhvpairs_a.size());

    for (auto la: lhvpairs_a) {
        m_accumuMargDistri[selected_vtx_a][la.first] += la.second;
    }

    for (auto lb: lhvpairs_b) {
        m_accumuMargDistri[selected_vtx_b][lb.first] += lb.second;
    }

    for (unsigned int i = 0; i < lhvpairs_a.size(); i++) {
        probarray[i] = lhvpairs_a[i].second;
    }
    m_accumuCondEntropy[selected_vtx_a] += entropy(probarray, (unsigned) (lhvpairs_a.size()));
    m_numAccumuCondEntropy[selected_vtx_a]++;

    for (unsigned int i = 0; i < lhvpairs_b.size(); i++) {
        probarray[i] = lhvpairs_b[i].second;
    }
    m_accumuCondEntropy[selected_vtx_b] += entropy(probarray, (unsigned) (lhvpairs_b.size()));
    m_numAccumuCondEntropy[selected_vtx_b]++;

}

Learner::~Learner() = default;

MutualInfo::~MutualInfo()
{
    for (unsigned int i=0; i < N_; ++i) {
        m_accumuMargDistri[i].clear();
    }
    m_accumuMargDistri.clear();

    m_accumuCondEntropy.clear();
    m_numAccumuCondEntropy.clear();
}
