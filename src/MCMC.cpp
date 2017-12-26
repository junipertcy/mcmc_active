#include "MCMC.h"
#include "Utility.cpp"


double MCMC::PI = 3.141592653589793;

MCMC::MCMC(TypeModel &tm, int blockmodeltype, bool groupcorrected) : m_typeModel(tm), m_blockModelType(blockmodeltype),
                                                                     m_groupCorrected(
                                                                             groupcorrected)//,m_bestTypeModel(tm)
{
    N_ = tm.getGraph().getNumVtx();
    Q_ = tm.getNumType();

    m_transProbSelect.resize(tm.getNumActiveType());

    m_maxsizeLogtable = N_ * N_ + 1;
    initLogTable(m_maxsizeLogtable);
    initLogGammaTable(2 + N_ * N_);

    dvtxClassifiMatrix = new double *[N_];
    for (unsigned int i = 0; i < N_; i++) {
        dvtxClassifiMatrix[i] = new double[Q_];
    }

    m_bestEdgeConnMatrix.resize(Q_);
    for (unsigned int i = 0; i < Q_; i++) {
        m_bestEdgeConnMatrix[i].resize(Q_);
    }
    m_bestGroupCardi.resize(Q_);
    m_bestVtxTypeTable.resize(N_);
    MAXLOGDOUBLE = log(numeric_limits<double>::max()) - 50;//-50 prevents the accumulation exceeding double limit
    MINLOGDOUBLE = log(numeric_limits<double>::min()) + 50;//

    m_LLHVariTable.resize(Q_);
}

unsigned MCMC::getTargetType(unsigned int mutateVtxNo) noexcept {
    calcLHVari(mutateVtxNo, m_typeModel);
    return calcTargetType();
}

void MCMC::calcLHVari(unsigned int vtxNo, TypeModel &typeModel) noexcept {
    m_LHVariPairs.clear();
    if (m_blockModelType == 1) {
        calcLHVariUDM1(vtxNo, typeModel);
    } else {
        cerr << "unrecognized model type! -- MCMC::calcLHVari()" << endl;
    }
}

//case 2.1: undirected for model type 1
void MCMC::calcLHVariUDM1(unsigned int v, TypeModel &typeModel) noexcept {
    unsigned int s, t, o;
    s = typeModel.getVtxType(v);
    unsigned int ns = typeModel.m_groupCardiTable[s];
    unsigned int nvs = typeModel.m_numTargetVtxGroup[v][s];

    double LHVari;
    pair<unsigned int, double> lhvariPair(s, 0.0);
    m_LHVariPairs.push_back(lhvariPair);
    m_LLHVariTable[s] = 0.0;
    for (t = 0; t < typeModel.getNumActiveType(); t++) {
        if (t == s)
            continue;
        unsigned int nt = typeModel.m_groupCardiTable[t];
        if (m_groupCorrected) {
            LHVari = (ns - 1) * getLog(ns - 1) + (nt + 1) * getLog(nt + 1) - ns * getLog(ns) -
                     nt * getLog(nt);//ns'*log(ns')+nt'*log(nt')-ns*log(ns)-nt*log(nt)
        } else {
            LHVari = 0.0;
        }
        for (o = 0; o < Q_; o++) {
            if (o == s || o == t)
                continue;
            //s<->o
            unsigned int eso = typeModel.m_numEdgesOf2Groups[s][o];
            unsigned int nvo = typeModel.m_numTargetVtxGroup[v][o];
            unsigned int no = typeModel.m_groupCardiTable[o];
            LHVari += getLogDivFac(eso - nvo, eso);
            LHVari += getLogDivFac(ns * no + 1, ns * no - no + 1);
            LHVari += getLogDivFac(ns * no - eso - no + nvo, ns * no - eso);
            //t<->o
            unsigned int eto = typeModel.m_numEdgesOf2Groups[t][o];
            LHVari += getLogDivFac(eto + nvo, eto);
            LHVari += getLogDivFac(nt * no + 1, nt * no + no + 1);
            LHVari += getLogDivFac(nt * no - eto + no - nvo, nt * no - eto);
        }
        //s<->t
        unsigned int est = typeModel.m_numEdgesOf2Groups[s][t];
        unsigned int nvt = typeModel.m_numTargetVtxGroup[v][t];

        /* Self-loop is not allowed */

        LHVari += getLogDivFac(est - nvt + nvs, est);
        LHVari += getLogDivFac(ns * nt + 1, (ns - 1) * (nt + 1) + 1);
        LHVari += getLogDivFac((ns - 1) * (nt + 1) - est + nvt - nvs, ns * nt - est);
        //s<->s
        unsigned int ess = typeModel.m_numEdgesOf2Groups[s][s];
        LHVari += getLogDivFac(ess - nvs, ess);
        LHVari += getLogDivFac(ns * (ns - 1) / 2 + 1, (ns - 1) * (ns - 2) / 2 + 1);
        LHVari += getLogDivFac((ns - 1) * (ns - 2) / 2 - ess + nvs, ns * (ns - 1) / 2 - ess);

        //t<->t
        unsigned int ett = typeModel.m_numEdgesOf2Groups[t][t];
        LHVari += getLogDivFac(ett + nvt, ett);
        LHVari += getLogDivFac(nt * (nt - 1) / 2 + 1, (nt + 1) * nt / 2 + 1);
        LHVari += getLogDivFac((nt + 1) * nt / 2 - ett - nvt, nt * (nt - 1) / 2 - ett);
        //
        pair<unsigned int, double> lhvariPair_(t, LHVari);
        m_LHVariPairs.push_back(lhvariPair_);
        m_LLHVariTable[t] = LHVari;
    }
}

unsigned MCMC::calcTargetType() noexcept {
    unsigned int i;
    double dsum = 0.0;
    double maxLogDif = m_LHVariPairs[0].second;
    double minLogDif = m_LHVariPairs[0].second;
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
    for (i = 0; i < m_LHVariPairs.size(); i++)
        m_LHVariPairs[i].second /= dsum;
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
double MCMC::getLogDivFac(unsigned int a, unsigned int b) {
    if (a > m_maxsizeLogtable || b > m_maxsizeLogtable) {
        cerr << "need more memory for logtable--getLogDivFac()" << endl;
//        delete[] m_logfactable;  // TODO: how to delete float_vec_t properly? Is it even needed?
        m_maxsizeLogtable = max(a, b) + 1000;
        initLogTable(m_maxsizeLogtable);
    }
    return m_logfactable[a] - m_logfactable[b];
}

// return  ln(a!)
double MCMC::getLogFac(unsigned int a) {
    if (a > m_maxsizeLogtable) {
        cerr << "need more memory for logtable -- getLogFac()\t" << a << endl;
//        delete[] m_logfactable;  //TODO: proper delete??
        m_maxsizeLogtable = a + 1000;
        initLogTable(m_maxsizeLogtable);
    }
    return m_logfactable[a];
}

double MCMC::getLog(unsigned int a) {
    if (a > m_maxsizeLogtable) {
        cerr << "need more memory for logtable -- getLog()\t" << a << endl;
//        delete[] m_logtable; //TODO
        m_maxsizeLogtable = a + 1000;
        initLogTable(m_maxsizeLogtable);
    }
    return m_logtable[a];
}

void MCMC::initVtxClassifiMatrix() {
    for (unsigned int i = 0; i < N_; ++i) {
        for (unsigned int j = 0; j < Q_; ++j) {
            dvtxClassifiMatrix[i][j] = 0.0;
        }
    }
}

void MCMC::updateVtxClassifiMatrix() {//should be called after mutateTypeModel()
    for (auto const &i: m_LHVariPairs) {
        dvtxClassifiMatrix[m_mutateVtxNo][i.first] += i.second;
    }
}

double **MCMC::getVtxClassifiMatrix() {
    double dSum;
    for (unsigned int i = 0; i < N_; i++) {
        dSum = 0.0;
        for (unsigned int j = 0; j < Q_; j++) {
            dSum += dvtxClassifiMatrix[i][j];
        }
        if (dSum != 0.0) {
            for (unsigned int j = 0; j < Q_; j++)
                dvtxClassifiMatrix[i][j] /= dSum;
        }
    }
    return dvtxClassifiMatrix;
}

void MCMC::randInitTypeModel(const set<unsigned> &topVtxSet) {
    m_typeModel.randInitGroups(topVtxSet);
    initLLHValue();
}

void MCMC::initBestTypeModel() {
    unsigned i, j;
    m_bestLLHvalue = calcLikelihood(m_typeModel);
    for (i = 0; i < Q_; i++) {
        for (j = 0; j < Q_; j++) {
            m_bestEdgeConnMatrix[i][j] = m_typeModel.m_numEdgesOf2Groups[i][j];
        }
        m_bestGroupCardi[i] = m_typeModel.m_groupCardiTable[i];
    }
    for (i = 0; i < N_; i++) {
        m_bestVtxTypeTable[i] = m_typeModel.getVtxType(i);
    }
}

void MCMC::updateBestTypeModel() {
    if (m_logLikelihoodValue > m_bestLLHvalue) {
        m_bestLLHvalue = calcLikelihood(m_typeModel);
        for (unsigned int i = 0; i < Q_; i++) {
            for (unsigned int j = 0; j < Q_; j++) {
                m_bestEdgeConnMatrix[i][j] = m_typeModel.m_numEdgesOf2Groups[i][j];
            }
            m_bestGroupCardi[i] = m_typeModel.m_groupCardiTable[i];
        }
        for (unsigned int i = 0; i < N_; i++) {
            m_bestVtxTypeTable[i] = m_typeModel.getVtxType(i);
        }
    }
}

void MCMC::mutateTypeModel(unsigned int v, unsigned int t) {
    m_mutateVtxNo = v;
    m_targetType = t;
    m_sourceType = m_typeModel.getVtxType(m_mutateVtxNo);
    if (m_targetType != m_sourceType) {
        m_typeModel.mutate(v, t);
        updateLLHValue();
        if (rand01() < 0.001) {
            double llhvalue = calcLikelihood(m_typeModel);
            if (fabs(llhvalue - m_logLikelihoodValue) > 0.00001) {
                cout << "likelihood calculation error! -- MCMC::mutateTypeModel()" << endl;
                cout << llhvalue << "\t" << m_logLikelihoodValue << endl;
                getchar();
            } else {
                m_logLikelihoodValue = llhvalue;
            }
        }
    }
}

double MCMC::calcLikelihood(const TypeModel &typemodel) {
    double llhvalue;
    if (m_blockModelType == 1) {
        llhvalue = calcLikelihoodM1(typemodel);
    } else {
        cerr << "unrecognized model type!" << endl;
        return -1;
    }
    return llhvalue;//-m_initllhvalue;
}

//This method returns the likelihood of the current Type Model (model type is 1)
double MCMC::calcLikelihoodM1(const TypeModel &typemodel) {
    unsigned int i, j, a, b;
    unsigned int numtype = typemodel.getNumType();
    double log_likelihood_value = 0.0;

    // Only for non-directed and no self-loop
    for (i = 0; i < numtype; i++) {
        for (j = i; j < numtype; j++) {
            a = typemodel.m_numEdgesOf2Groups[i][j];
            if (i == j)
                b = typemodel.m_groupCardiTable[i] * (typemodel.m_groupCardiTable[i] - 1) / 2 - a;
            else
                b = typemodel.m_groupCardiTable[i] * typemodel.m_groupCardiTable[j] - a;
            log_likelihood_value = log_likelihood_value + getLogFac(a) + getLogFac(b) - getLogFac(a + b + 1);
        }
    }
    return log_likelihood_value;
}

