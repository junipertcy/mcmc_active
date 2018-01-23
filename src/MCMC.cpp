#include "MCMC.h"
#include "Utility.cpp"
#include "output_functions.h"


double MCMC::PI = 3.141592653589793;

MCMC::MCMC(TypeModel &tm, unsigned int blockmodeltype, bool groupcorrected):
        m_typeModel(tm),
        m_blockModelType(blockmodeltype),
        m_groupCorrected(groupcorrected)
{
    N_ = tm.getGraph().get_N();
    Q_ = tm.get_Q();

    m_transProbSelect.resize(tm.getNumActiveType());

    m_maxsizeLogtable = N_ * N_ + 1;
    initLogTable(m_maxsizeLogtable);
    initLogGammaTable(2 + N_ * N_);

    dvtxClassifiMatrix = new double *[N_];
    for (unsigned int i = 0; i < N_; i++) {
        dvtxClassifiMatrix[i] = new double[Q_];
    }

    MAXLOGDOUBLE = log(numeric_limits<double>::max()) - 50;//-50 prevents the accumulation exceeding double limit
    MINLOGDOUBLE = log(numeric_limits<double>::min()) + 50;//

    m_LLHVariTable.resize(Q_, 0.);
}

unsigned int MCMC::get_selected_vtx() const noexcept { return m_mutateVtxNo; }

const TypeModel &MCMC::getTypeModel() const noexcept { return m_typeModel; }

const vector<pair<unsigned int,double>> &MCMC::get_likelihood_variation_pairs() const noexcept { return m_LHVariPairs; }

unsigned int MCMC::getTargetType(unsigned int mutateVtxNo) noexcept {
    compute_likelihood(mutateVtxNo, m_typeModel);
    return compute_gibbs_jump();
}

void MCMC::compute_likelihood(unsigned int vtxNo, TypeModel &typeModel) noexcept {
    m_LHVariPairs.clear();
    if (m_blockModelType == 1) {
        _compute_likelihood(vtxNo, typeModel);
    } else {
        cerr << "unrecognized model type! -- MCMC::calcLHVari()" << endl;
    }
}

//case 2.1: undirected for model type 1
void MCMC::_compute_likelihood(unsigned int vtx, TypeModel &typeModel) noexcept {
    unsigned int source_group, target_group, q;
    source_group = typeModel.get_membership(vtx);
    unsigned int ns = typeModel.n_r_[source_group];
    unsigned int nvs = typeModel.m_numTargetVtxGroup[vtx][source_group];

    double LHVari;
    pair<unsigned int, double> lhvariPair(source_group, 0.0);
    m_LHVariPairs.push_back(lhvariPair);

    m_LLHVariTable[source_group] = 0.0;
    for (target_group = 0; target_group < typeModel.getNumActiveType(); target_group++) {
        if (target_group == source_group) {
            continue;
        }
        unsigned int nt = typeModel.n_r_[target_group];
        if (m_groupCorrected) {
            LHVari = (ns - 1) * getLog(ns - 1) + (nt + 1) * getLog(nt + 1) - ns * getLog(ns) -
                     nt * getLog(nt);//ns'*log(ns')+nt'*log(nt')-ns*log(ns)-nt*log(nt)
        } else {
            LHVari = 0.0;
        }
        for (q = 0; q < Q_; q++) {
            if (q == source_group || q == target_group)
                continue;

            //source_group<->q
            unsigned int eso = typeModel.e_rs_[source_group][q];
            unsigned int nvo = typeModel.m_numTargetVtxGroup[vtx][q];
            unsigned int no = typeModel.n_r_[q];
            LHVari += getLogDivFac(eso - nvo, eso);
            LHVari += getLogDivFac(ns * no + 1, ns * no - no + 1);
            LHVari += getLogDivFac(ns * no - eso - no + nvo, ns * no - eso);

            //target_group<->q
            unsigned int eto = typeModel.e_rs_[target_group][q];
            LHVari += getLogDivFac(eto + nvo, eto);
            LHVari += getLogDivFac(nt * no + 1, nt * no + no + 1);
            LHVari += getLogDivFac(nt * no - eto + no - nvo, nt * no - eto);
        }
        //source_group<->target_group
        unsigned int est = typeModel.e_rs_[source_group][target_group];
        unsigned int nvt = typeModel.m_numTargetVtxGroup[vtx][target_group];

        /* Self-loop is not allowed */
        LHVari += getLogDivFac(est - nvt + nvs, est);
        LHVari += getLogDivFac(ns * nt + 1, (ns - 1) * (nt + 1) + 1);
        LHVari += getLogDivFac((ns - 1) * (nt + 1) - est + nvt - nvs, ns * nt - est);

        //source_group<->source_group
        unsigned int ess = typeModel.e_rs_[source_group][source_group];
        LHVari += getLogDivFac(ess - nvs, ess);
        LHVari += getLogDivFac(ns * (ns - 1) / 2 + 1, (ns - 1) * (ns - 2) / 2 + 1);
        LHVari += getLogDivFac((ns - 1) * (ns - 2) / 2 - ess + nvs, ns * (ns - 1) / 2 - ess);

        //target_group<->target_group
        unsigned int ett = typeModel.e_rs_[target_group][target_group];
        LHVari += getLogDivFac(ett + nvt, ett);
        LHVari += getLogDivFac(nt * (nt - 1) / 2 + 1, (nt + 1) * nt / 2 + 1);
        LHVari += getLogDivFac((nt + 1) * nt / 2 - ett - nvt, nt * (nt - 1) / 2 - ett);
        //
        pair<unsigned int, double> lhvariPair_(target_group, LHVari);
        m_LHVariPairs.push_back(lhvariPair_);
        m_LLHVariTable[target_group] = LHVari;
//        std::clog << "m_LHVariPairs.size() = " << m_LHVariPairs.size() << "\n";
//        std::clog << "source_group = " << m_LHVariPairs[0].first << "; LHVari = " << m_LHVariPairs[0].second <<"\n";
//        std::clog << "target_group = " << m_LHVariPairs[1].first << "; LHVari = " << m_LHVariPairs[1].second <<"\n";
//        std::clog << "---\n";
    }
}

unsigned int MCMC::compute_gibbs_jump() noexcept {
    unsigned int i;
    double dsum = 0.0;
    double maxLogDif = 0;
    double minLogDif = 0;

    for (i = 1; i < m_LHVariPairs.size(); i++) {
        if (m_LHVariPairs[i].second > maxLogDif)
            maxLogDif = m_LHVariPairs[i].second;
        if (m_LHVariPairs[i].second < minLogDif)
            minLogDif = m_LHVariPairs[i].second;
    }
    double newMaxLogDif = (maxLogDif - minLogDif) / 2;
    double offset = newMaxLogDif - maxLogDif;
    for (i = 0; i < m_LHVariPairs.size(); i++) {
        m_LHVariPairs[i].second += offset;
    }
    if (newMaxLogDif > MAXLOGDOUBLE) {
        double shrinkRatio = newMaxLogDif / MAXLOGDOUBLE;
        for (i = 0; i < m_LHVariPairs.size(); i++) {
            m_LHVariPairs[i].second /= shrinkRatio;
        }
    }
    for (i = 0; i < m_LHVariPairs.size(); i++) {
        m_LHVariPairs[i].second = exp(m_LHVariPairs[i].second);
        dsum += m_LHVariPairs[i].second;
        m_transProbSelect[i] = dsum;
    }
    for (i = 0; i < m_LHVariPairs.size(); i++) {
        m_LHVariPairs[i].second /= dsum;
    }

    // Choose which label to change to, for vertex vtx; heat-bath jump.
    int index = getIndexProb(m_transProbSelect, (int) m_LHVariPairs.size());
    unsigned int targetType = m_LHVariPairs[index].first;
    return targetType;
}

//m_logfactable[i]=ln(i!)
void MCMC::initLogTable(unsigned int size) noexcept {
    m_logfactable.resize(size + 1);
    m_logtable.resize(size + 1);
    m_logfactable[0] = 0;
    m_logtable[0] = 0;
    for (unsigned i = 1; i <= size; i++) {
        m_logfactable[i] = m_logfactable[i - 1] + log((double) i);
        m_logtable[i] = log((double) i);
    }
}

//log_gamma(0),log_gamma(0.5),log_gamma(1)...; log_gamma(i)=m_loggammatable[2i];
void MCMC::initLogGammaTable(unsigned int size) noexcept {
    m_loggammatable.resize(size + 1);
    m_loggammatable[0] = 0;  //dummy, will never be used;
    m_loggammatable[1] = log(sqrt(PI));  //log_gamma(0.5)
    m_loggammatable[2] = log(1.0);  //log_gamma(1)
    for (unsigned int i = 3; i < size; i++)
        m_loggammatable[i] = log((0.5 * i - 1)) + m_loggammatable[i - 2];//log_gamma(0.5i)
}

// return  ln(a!/b!)
double MCMC::getLogDivFac(unsigned int a, unsigned int b) noexcept {
    if (a > m_maxsizeLogtable || b > m_maxsizeLogtable) {
        cerr << "need more memory for logtable--getLogDivFac()" << endl;
//        delete[] m_logfactable;  // TODO: how to delete float_vec_t properly? Is it even needed?
        m_maxsizeLogtable = max(a, b) + 1000;
        initLogTable(m_maxsizeLogtable);
    }
    return m_logfactable[a] - m_logfactable[b];
}

// return  ln(a!)
double MCMC::getLogFac(unsigned int a) noexcept {
    if (a > m_maxsizeLogtable) {
        cerr << "need more memory for logtable -- getLogFac()\t" << a << endl;
//        delete[] m_logfactable;  //TODO: proper delete??
        m_maxsizeLogtable = a + 1000;
        initLogTable(m_maxsizeLogtable);
    }
    return m_logfactable[a];
}

double MCMC::getLog(unsigned int a) noexcept {
    if (a > m_maxsizeLogtable) {
        cerr << "need more memory for logtable -- getLog()\t" << a << endl;
//        delete[] m_logtable; //TODO
        m_maxsizeLogtable = a + 1000;
        initLogTable(m_maxsizeLogtable);
    }
    return m_logtable[a];
}

void MCMC::randInitTypeModel(const set<unsigned> &topVtxSet) noexcept {
    m_typeModel.randInitGroups(topVtxSet);
    init_log_likelihood();
}


void MCMC::gibbs_jump(unsigned int vtx, unsigned int t) noexcept {
    m_mutateVtxNo = vtx;
    m_targetType = t;
    m_sourceType = m_typeModel.get_membership(m_mutateVtxNo);
    if (m_targetType != m_sourceType) {
        m_typeModel.mutate(vtx, t);

        update_log_likelihood();

        if (rand01() < 0.001) {
            double llhvalue = calcLikelihood(m_typeModel);
            if (fabs(llhvalue - m_logLikelihoodValue) > 0.00001) {
                cout << "likelihood calculation error! -- MCMC::gibbs_jump()" << endl;
                cout << llhvalue << "\t" << m_logLikelihoodValue << endl;
                getchar();
            } else {
                m_logLikelihoodValue = llhvalue;
            }
        }
    }
}

double MCMC::calcLikelihood(const TypeModel &typemodel) noexcept {
    double llhvalue;
    if (m_blockModelType == 1) {
        llhvalue = calcLikelihoodM1(typemodel);
    } else {
        cerr << "unrecognized model type!" << endl;
        return -1;
    }
    return llhvalue;  //-m_initllhvalue;
}

//This method returns the likelihood of the current Type Model (model type is 1)
double MCMC::calcLikelihoodM1(const TypeModel &typemodel) noexcept {
    double log_likelihood_value = 0.0;

    // Only for non-directed and no self-loop
    for (unsigned int i = 0; i < Q_; i++) {
        for (unsigned int j = i; j < Q_; j++) {
            unsigned int a = typemodel.e_rs_[i][j];
            unsigned int b;
            if (i == j) {
                b = typemodel.n_r_[i] * (typemodel.n_r_[i] - 1) / 2 - a;
            } else {
                b = typemodel.n_r_[i] * typemodel.n_r_[j] - a;
            }
            log_likelihood_value = log_likelihood_value + getLogFac(a) + getLogFac(b) - getLogFac(a + b + 1);
        }
    }
    return log_likelihood_value;
}

void MCMC::init_log_likelihood() noexcept { m_logLikelihoodValue = calcLikelihood(m_typeModel); }

void MCMC::update_log_likelihood() noexcept {
    m_logLikelihoodValue += m_LLHVariTable[m_targetType];
}
