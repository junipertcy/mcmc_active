#include "MCMC.h"
#include "Utility.cpp"


double MCMC::PI = 3.141592653589793;

MCMC::MCMC(TypeModel &tm, int blockmodeltype, bool groupcorrected) : m_typeModel(tm), m_blockModelType(blockmodeltype),
                                                                     m_groupCorrected(
                                                                             groupcorrected)//,m_bestTypeModel(tm)
{
    unsigned int i;
    m_numVtx = tm.getGraph().getNumVtx();
    m_numType = tm.getNumType();

    m_transProbSelect.resize(tm.getNumActiveType());

    m_maxsizeLogtable = m_numVtx * m_numVtx + 1;
    initLogTable(m_maxsizeLogtable);
    initLogGammaTable(2 + m_numVtx * m_numVtx);

    dvtxClassifiMatrix = new double *[m_numVtx];
    for (i = 0; i < m_numVtx; i++) {
        dvtxClassifiMatrix[i] = new double[m_numType];
    }

    m_bestEdgeConnMatrix.resize(m_numType);
    for (i = 0; i < m_numType; i++) {
        m_bestEdgeConnMatrix[i].resize(m_numType);
    }
    m_bestGroupCardi.resize(m_numType);
    m_bestVtxTypeTable.resize(m_numVtx);
    MAXLOGDOUBLE = log(numeric_limits<double>::max()) - 50;//-50 prevents the accumulation exceeding double limit
    MINLOGDOUBLE = log(numeric_limits<double>::min()) + 50;//

    m_LLHVariTable.resize(m_numType);
}

unsigned MCMC::getTargetType(unsigned int mutateVtxNo) noexcept {
    calcLHVari(mutateVtxNo, m_typeModel);
    return calcTargetType();
}

void MCMC::calcLHVari(unsigned int vtxNo, TypeModel &typeModel) noexcept {

    m_LHVariPairs.clear();
    bool isDirected = typeModel.getGraph().isDirected();
//    std::clog << "======== isDirected is : " << !isDirected << "\n";
    if (!isDirected && m_blockModelType == 1) {
        calcLHVariUDM1(vtxNo, typeModel);
    } else {
        cerr << "unrecognized model type! -- MCMC::calcLHVari()" << endl;
    }
}

//case 2.1: undirected for model type 1
void MCMC::calcLHVariUDM1(unsigned int v, TypeModel &typeModel) noexcept {
    unsigned s, t, o;
    s = typeModel.getVtxType(v);
    unsigned ns = typeModel.m_groupCardiTable[s];
    unsigned nvs = typeModel.m_numTargetVtxGroup[v][s];
    bool hasSelfLoop = typeModel.getGraph().hasSelfloop();
    unsigned lvv = 0;
    if (typeModel.getGraph().vtxHasSelfloop(v))
        lvv = 1;
    double LHVari;
    pair<unsigned, double> lhvariPair(s, 0.0);
    m_LHVariPairs.push_back(lhvariPair);
    m_LLHVariTable[s] = 0.0;
    for (t = 0; t < typeModel.getNumActiveType(); t++) {
        if (t == s)
            continue;
        unsigned nt = typeModel.m_groupCardiTable[t];
        if (m_groupCorrected)
            LHVari = (ns - 1) * getLog(ns - 1) + (nt + 1) * getLog(nt + 1) - ns * getLog(ns) -
                     nt * getLog(nt);//ns'*log(ns')+nt'*log(nt')-ns*log(ns)-nt*log(nt)
        else
            LHVari = 0.0;
        for (o = 0; o < m_numType; o++) {
            if (o == s || o == t)
                continue;
            //s<->o
            unsigned eso = typeModel.m_numEdgesOf2Groups[s][o];
            unsigned nvo = typeModel.m_numTargetVtxGroup[v][o];
            unsigned no = typeModel.m_groupCardiTable[o];
            LHVari += getLogDivFac(eso - nvo, eso);
            LHVari += getLogDivFac(ns * no + 1, ns * no - no + 1);
            LHVari += getLogDivFac(ns * no - eso - no + nvo, ns * no - eso);
            //t<->o
            unsigned eto = typeModel.m_numEdgesOf2Groups[t][o];
            LHVari += getLogDivFac(eto + nvo, eto);
            LHVari += getLogDivFac(nt * no + 1, nt * no + no + 1);
            LHVari += getLogDivFac(nt * no - eto + no - nvo, nt * no - eto);
        }
        //s<->t
        unsigned est = typeModel.m_numEdgesOf2Groups[s][t];
        unsigned nvt = typeModel.m_numTargetVtxGroup[v][t];
        if (hasSelfLoop) {
            LHVari += getLogDivFac(est - nvt + nvs - lvv, est);
            LHVari += getLogDivFac(ns * nt + 1, (ns - 1) * (nt + 1) + 1);
            LHVari += getLogDivFac((ns - 1) * (nt + 1) - est + nvt - nvs + lvv, ns * nt - est);
        } else {
            LHVari += getLogDivFac(est - nvt + nvs, est);
            LHVari += getLogDivFac(ns * nt + 1, (ns - 1) * (nt + 1) + 1);
            LHVari += getLogDivFac((ns - 1) * (nt + 1) - est + nvt - nvs, ns * nt - est);
        }
        //s<->s
        unsigned ess = typeModel.m_numEdgesOf2Groups[s][s];
        if (hasSelfLoop) {
            LHVari += getLogDivFac(ess - nvs, ess);
            LHVari += getLogDivFac(ns * (ns + 1) / 2 + 1, ns * (ns - 1) / 2 + 1);
            LHVari += getLogDivFac(ns * (ns - 1) / 2 - ess + nvs, ns * (ns + 1) / 2 - ess);
        } else {
            LHVari += getLogDivFac(ess - nvs, ess);
            LHVari += getLogDivFac(ns * (ns - 1) / 2 + 1, (ns - 1) * (ns - 2) / 2 + 1);
            LHVari += getLogDivFac((ns - 1) * (ns - 2) / 2 - ess + nvs, ns * (ns - 1) / 2 - ess);
        }
        //t<->t
        unsigned ett = typeModel.m_numEdgesOf2Groups[t][t];
        if (hasSelfLoop) {
            LHVari += getLogDivFac(ett + nvt + lvv, ett);
            LHVari += getLogDivFac(nt * (nt + 1) / 2 + 1, (nt + 1) * (nt + 2) / 2 + 1);
            LHVari += getLogDivFac((nt + 1) * (nt + 2) / 2 - ett - nvt - lvv, nt * (nt + 1) / 2 - ett);
        } else {
            LHVari += getLogDivFac(ett + nvt, ett);
            LHVari += getLogDivFac(nt * (nt - 1) / 2 + 1, (nt + 1) * nt / 2 + 1);
            LHVari += getLogDivFac((nt + 1) * nt / 2 - ett - nvt, nt * (nt - 1) / 2 - ett);
        }
        //
        pair<unsigned, double> lhvariPair(t, LHVari);
        m_LHVariPairs.push_back(lhvariPair);
        m_LLHVariTable[t] = LHVari;
    }
}

unsigned MCMC::calcTargetType() noexcept {
    unsigned i;
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
    unsigned index = getIndexProb(m_transProbSelect, m_LHVariPairs.size());
    unsigned targetType = m_LHVariPairs[index].first;
    return targetType;
}

//m_logfactable[i]=ln(i!)
void MCMC::initLogTable(unsigned size) noexcept {
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
void MCMC::initLogGammaTable(unsigned size) noexcept {
    unsigned i;
    m_loggammatable.resize(size + 1);
    m_loggammatable[0] = 0;//dummy, will never be used;
    m_loggammatable[1] = log(sqrt(PI));//log_gamma(0.5)
    m_loggammatable[2] = log(1.0);//log_gamma(1)
    for (i = 3; i < size; i++)
        m_loggammatable[i] = log((0.5 * i - 1)) + m_loggammatable[i - 2];//log_gamma(0.5i)
}

// return  ln(a!/b!)
double MCMC::getLogDivFac(unsigned a, unsigned b) {
    if (a > m_maxsizeLogtable || b > m_maxsizeLogtable) {
        cerr << "need more memory for logtable--getLogDivFac()" << endl;
//        delete[] m_logfactable;  // TODO: how to delete float_vec_t properly? Is it even needed?
        m_maxsizeLogtable = max(a, b) + 1000;
        initLogTable(m_maxsizeLogtable);
    }
    return m_logfactable[a] - m_logfactable[b];
}

// return  ln(a!)
double MCMC::getLogFac(unsigned a) {
    if (a > m_maxsizeLogtable) {
        cerr << "need more memory for logtable -- getLogFac()\t" << a << endl;
//        delete[] m_logfactable;  //TODO: proper delete??
        m_maxsizeLogtable = a + 1000;
        initLogTable(m_maxsizeLogtable);
    }
    return m_logfactable[a];
}

double MCMC::getLog(unsigned a) {
    if (a > m_maxsizeLogtable) {
        cerr << "need more memory for logtable -- getLog()\t" << a << endl;
//        delete[] m_logtable; //TODO
        m_maxsizeLogtable = a + 1000;
        initLogTable(m_maxsizeLogtable);
    }
    return m_logtable[a];
}

void MCMC::initVtxClassifiMatrix() {
    unsigned i, j;
    for (i = 0; i < m_numVtx; i++) {
        for (j = 0; j < m_numType; j++) {
            dvtxClassifiMatrix[i][j] = 0.0;
        }
    }
}

void MCMC::updateVtxClassifiMatrix() {//should be called after mutateTypeModel()
    unsigned i;
    for (i = 0; i < m_LHVariPairs.size(); i++)
        dvtxClassifiMatrix[m_mutateVtxNo][m_LHVariPairs[i].first] += m_LHVariPairs[i].second;
}

double **MCMC::getVtxClassifiMatrix() {
    unsigned i, j;
    double dSum;
    for (i = 0; i < m_numVtx; i++) {
        dSum = 0.0;
        for (j = 0; j < m_numType; j++) {
            dSum += dvtxClassifiMatrix[i][j];
        }
        if (dSum != 0) {
            for (j = 0; j < m_numType; j++)
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
    //m_bestTypeModel=m_typeModel;
    for (i = 0; i < m_numType; i++) {
        for (j = 0; j < m_numType; j++) {
            m_bestEdgeConnMatrix[i][j] = m_typeModel.m_numEdgesOf2Groups[i][j];
        }
        m_bestGroupCardi[i] = m_typeModel.m_groupCardiTable[i];
    }
    for (i = 0; i < m_numVtx; i++) {
        m_bestVtxTypeTable[i] = m_typeModel.getVtxType(i);
    }
}

void MCMC::updateBestTypeModel() {
    unsigned i, j;
    if (m_logLikelihoodValue > m_bestLLHvalue) {
        m_bestLLHvalue = calcLikelihood(m_typeModel);
        for (i = 0; i < m_numType; i++) {
            for (j = 0; j < m_numType; j++) {
                m_bestEdgeConnMatrix[i][j] = m_typeModel.m_numEdgesOf2Groups[i][j];
            }
            m_bestGroupCardi[i] = m_typeModel.m_groupCardiTable[i];
        }
        for (i = 0; i < m_numVtx; i++) {
            m_bestVtxTypeTable[i] = m_typeModel.getVtxType(i);
        }
    }
}

void MCMC::mutateTypeModel(unsigned v, unsigned t) {
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
    unsigned i, j, a, b;
    unsigned numtype = typemodel.getNumType();
    //double likelihood_value=0.0;
    double log_likelihood_value = 0.0;
    if (typemodel.getGraph().isDirected() && typemodel.getGraph().hasSelfloop()) {
        for (i = 0; i < numtype; i++) {
            for (j = 0; j < numtype; j++) {
                a = typemodel.m_numEdgesOf2Groups[i][j];
                b = typemodel.m_groupCardiTable[i] * typemodel.m_groupCardiTable[j] - a;
                log_likelihood_value = log_likelihood_value + getLogFac(a) + getLogFac(b) - getLogFac(a + b + 1);
            }
        }
    }
    if (typemodel.getGraph().isDirected() && !typemodel.getGraph().hasSelfloop()) {
        for (i = 0; i < numtype; i++) {
            for (j = 0; j < numtype; j++) {
                a = typemodel.m_numEdgesOf2Groups[i][j];
                if (i == j)
                    b = typemodel.m_groupCardiTable[i] * (typemodel.m_groupCardiTable[i] - 1) - a;
                else
                    b = typemodel.m_groupCardiTable[i] * typemodel.m_groupCardiTable[j] - a;
                log_likelihood_value = log_likelihood_value + getLogFac(a) + getLogFac(b) - getLogFac(a + b + 1);
            }
        }
    }
    if (!typemodel.getGraph().isDirected() && typemodel.getGraph().hasSelfloop()) {
        for (i = 0; i < numtype; i++) {
            for (j = i; j < numtype; j++) {
                a = typemodel.m_numEdgesOf2Groups[i][j];
                if (i == j)
                    b = typemodel.m_groupCardiTable[i] * (typemodel.m_groupCardiTable[i] - 1) / 2 +
                        typemodel.m_groupCardiTable[i] - a;
                else
                    b = typemodel.m_groupCardiTable[i] * typemodel.m_groupCardiTable[j] - a;
                log_likelihood_value = log_likelihood_value + getLogFac(a) + getLogFac(b) - getLogFac(a + b + 1);
            }
        }
    }
    if (!typemodel.getGraph().isDirected() && !typemodel.getGraph().hasSelfloop()) {
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
    }
    return log_likelihood_value;
}

