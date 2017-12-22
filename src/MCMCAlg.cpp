#include "MCMCAlg.h"
#include "Utility.cpp"

double MCMCAlg::stabRatio = 0.5;
unsigned MCMCAlg::numAccuracyBlock = 10;

MCMCAlg::MCMCAlg(const Graph &graph, unsigned numTypeInModel, set<unsigned> &frozentypes, unsigned numOptInit,
                 long numOptStep, unsigned numLearnerInit, long numLearnerStep, unsigned numPhase, unsigned numTop,
                 unsigned learningMethod, unsigned modelType, bool groupcorrected):
       m_graph(graph) {

    m_typeModelA = std::make_unique<TypeModel>(graph, numTypeInModel, frozentypes);
    m_typeModelB = std::make_unique<TypeModel>(graph, numTypeInModel, frozentypes);
    m_MC_A = std::make_unique<MCMC>(*m_typeModelA, modelType, groupcorrected);
    m_MC_B = std::make_unique<MCMC>(*m_typeModelB, modelType, groupcorrected);
    unsigned i;
    m_numOptInit = numOptInit;
    m_numOptStep = numOptStep;
    m_numLearnerInit = numLearnerInit;
    m_numLearnerStep = numLearnerStep;
    m_numPhase = numPhase;
    m_numTop = numTop;
    m_method = learningMethod;  // indicates which active learner is used: 1-MutualInfo, 2-AvgAgree, 3-RandomLearner, 4-TopSeqLearner, 5-MaxDegree, 6-MaxBtwn.

    switch (m_method) {
        default:
            m_learner = std::make_unique<MutualInfo>(*m_MC_A, *m_MC_B, m_topVtxSet);
    }
    m_alltypefixed = false;

    m_numVtx = m_typeModelA->getGraph().getNumVtx();
    m_numTypeGraph = m_typeModelA->getGraph().getNumType();
    m_numTypeModel = m_typeModelA->getNumType();
    for (i = 0; i < m_numVtx; i++) {
        m_remainVtxSet.insert(i);
    }
    //init top vtx seq
    for (i = 0; i < m_numVtx; i++) {
        unsigned gtype = graph.getVertex(i).getType();
        if (m_typeModelA->isGTypeFrozen(gtype)) {
            m_topVtxSet.insert(i);
            m_remainVtxSet.erase(i);
            m_topVtxSeq.push_back(i);
        }
    }
    //showiset(m_topVtxSet);
    m_bestClfcInPhase.resize(m_numVtx);
    m_accumuMargDistri.resize(m_numVtx);
    for (i = 0; i < m_numVtx; i++) {
        m_accumuMargDistri[i].resize(m_numTypeModel);
    }

    if (m_method == 3 || m_method == 4 || m_method == 5) {
        m_numLearnerInit = 0;
        m_numLearnerStep = 0;
    }

}

void MCMCAlg::runMCMCAlg() {
    initAccuracyMatrix();
    m_curPhase = 1;
#ifdef RESUME
    m_curPhase=m_resumePhase;
#endif
    for (; (m_curPhase <= m_numPhase) && (m_remainVtxSet.size() > 0); m_curPhase++) {
        runOnePhase();
    }
}

void MCMCAlg::runOnePhase() {

    unsigned i, j;

    checkTypeFixed();

    m_dSumLHvalueInPhase = 0.0;

    m_lNumLHvalueInPhase = 0;

    m_learner->resetForNextPhase();

    m_bestLLHvalueInPhase = LARGE_NEGATIVE_VALUE;

    initAccumuMargDistri();

    (this->m_MC_A)->initGroupConnMatrix();
    (this->m_MC_B)->initGroupConnMatrix();
    (this->m_MC_A)->initVtxClassifiMatrix();
    (this->m_MC_B)->initVtxClassifiMatrix();


    for (m_curInit = 1; m_curInit <= max(m_numOptInit, m_numLearnerInit); m_curInit++) {
        runOneInit();
        updateBestTypeModelInPhase();
    }

    updateAccuracyMatrix();

    getTopVtx();

    double **typeClassifiMatrix = new double *[m_numTypeGraph];
    for (i = 0; i < m_numTypeGraph; i++)
        typeClassifiMatrix[i] = new double[m_numTypeModel];
    for (i = 0; i < m_numTypeGraph; i++) {
        for (j = 0; j < m_numTypeModel; j++) {
            typeClassifiMatrix[i][j] = 0.0;
        }
    }
    int realtype;
    double **vtxClassifiMatrixA = m_MC_A->getVtxClassifiMatrix();
    double **vtxClassifiMatrixB = m_MC_B->getVtxClassifiMatrix();
    for (i = 0; i < m_numVtx; i++) {
        realtype = m_MC_A->getTypeModel().getGraph().getVertex(i).getType();
        for (j = 0; j < m_numTypeModel; j++) {
            typeClassifiMatrix[realtype][j] +=/*vtxClassifiMatrixA[i][j];//*/
                    (vtxClassifiMatrixA[i][j] + vtxClassifiMatrixB[i][j]) / 2;
        }
    }

    unsigned *topvtxGroupCardi = new unsigned[m_numTypeGraph];
    for (i = 0; i < m_numTypeGraph; i++)
        topvtxGroupCardi[i] = 0;
    set<unsigned>::const_iterator csiter;
    for (csiter = m_topVtxSet.begin(); csiter != m_topVtxSet.end(); csiter++) {
        realtype = (m_typeModelA->m_graph).getVertex(*csiter).getType();
        topvtxGroupCardi[realtype]++;
    }
    //
    for (i = 0; i < m_numTypeGraph; i++)
        delete[] typeClassifiMatrix[i];
    delete[] typeClassifiMatrix;
    delete[] topvtxGroupCardi;
}

void MCMCAlg::runOneInit() {


    m_MC_A->randInitTypeModel(this->m_topVtxSet);
    m_MC_B->randInitTypeModel(this->m_topVtxSet);

    (this->m_MC_A)->initBestTypeModel();
    (this->m_MC_B)->initBestTypeModel();

    m_learner->resetForNextInit();//this is very important for average agreement learner.

    long maxstep = max(m_numOptStep, m_numLearnerStep);
    if (min(m_numOptInit, m_numLearnerInit) < m_curInit) {
        if (m_numOptInit < m_numLearnerInit)
            maxstep = m_numLearnerStep;
        else
            maxstep = m_numOptStep;
    }

    for (m_curStep = 1; m_curStep <= maxstep; m_curStep++) {
        runOneStep();
    }
}

void MCMCAlg::runOneStep() {
    unsigned mutateVtxNo;
    unsigned sourceType;
    unsigned targetType;
    //typeModelA
    do {
        mutateVtxNo = getRandRemainVtx();
        sourceType = m_typeModelA->m_vtxTypeTable[mutateVtxNo];
    } while (m_typeModelA->m_groupCardiTable[sourceType] == 1);
    targetType = m_MC_A->getTargetType(mutateVtxNo);

    m_MC_A->mutateTypeModel(mutateVtxNo, targetType);


    if ((m_curInit <= m_numOptInit) && (m_curStep > stabRatio * double(m_numOptStep)) && (m_curStep <= m_numOptStep)) {
        updateAccumuMargDistri(mutateVtxNo, *m_MC_A);
        (this->m_MC_A)->updateGroupConnMatrix();
        (this->m_MC_A)->updateVtxClassifiMatrix();
        (this->m_MC_A)->updateBestTypeModel();
        m_dSumLHvalueInPhase += m_MC_A->getLHvalue();
        this->m_lNumLHvalueInPhase++;
    }

    //typeModelB
    do {
        mutateVtxNo = getRandRemainVtx();
        sourceType = m_typeModelB->m_vtxTypeTable[mutateVtxNo];
    } while (m_typeModelB->m_groupCardiTable[sourceType] == 1);
    targetType = m_MC_B->getTargetType(mutateVtxNo);

    m_MC_B->mutateTypeModel(mutateVtxNo, targetType);

    if ((m_curInit <= m_numOptInit) && (m_curStep > stabRatio * double(m_numOptStep)) && (m_curStep <= m_numOptStep)) {
        updateAccumuMargDistri(mutateVtxNo, *m_MC_B);
        (this->m_MC_B)->updateGroupConnMatrix();
        (this->m_MC_B)->updateVtxClassifiMatrix();
        (this->m_MC_B)->updateBestTypeModel();
        m_dSumLHvalueInPhase += m_MC_B->getLHvalue();
        this->m_lNumLHvalueInPhase++;
    }

    //
    if ((m_curInit <= m_numLearnerInit) && (m_curStep > stabRatio * double(m_numLearnerStep)) &&
        (m_curStep <= m_numLearnerStep)) {
        m_learner->updateData();
    }
}

int MCMCAlg::getRandRemainVtx() {
    unsigned index = randN(m_remainVtxSet.size());
    set<unsigned>::iterator setiter = m_remainVtxSet.begin();
    for (unsigned i = 0; i < index; i++)
        setiter++;
    return *setiter;
}

void MCMCAlg::getTopVtx() {
    unsigned i;
    unsigned *arrayTop = new unsigned[m_numVtx];//may output more than "m_numTop" vertices
    unsigned int numtop = m_learner->getTopVtx(arrayTop, m_numTop);
    for (i = 0; i < numtop; i++) {
        m_topVtxSet.insert(arrayTop[i]);
        m_remainVtxSet.erase(arrayTop[i]);
        m_topVtxSeq.push_back(arrayTop[i]);
    }
    delete[] arrayTop;
}

void MCMCAlg::initAccuracyMatrix() {
    m_accuracyMatrix.resize(numAccuracyBlock - 1);
    for (unsigned int i = 0; i < numAccuracyBlock - 1; i++)
        m_accuracyMatrix[i].resize(m_numPhase);
    for (unsigned int i = 0; i < numAccuracyBlock - 1; i++)
        for (unsigned int j = 0; j < m_numPhase; j++)
            m_accuracyMatrix[i][j] = 0.0;
}

void MCMCAlg::updateAccuracyMatrix() {
    unsigned int count = 0;
    double accuracyRatio = 1.0 / numAccuracyBlock;
    for (unsigned int vtxno = 0; vtxno < m_numVtx; vtxno++) {
        if (m_topVtxSet.count(vtxno)) {
            continue;
        }
        unsigned int realtype = m_typeModelA->getGraph().getVertex(vtxno).getType();
        double accuracyvalue = m_accumuMargDistri[vtxno][realtype];
        for (unsigned int i = 1; i < numAccuracyBlock; i++) {
            if (accuracyvalue >= accuracyRatio * i) {
                m_accuracyMatrix[i - 1][m_curPhase - 1]++;
            } else break;
        }
        count++;
    }
    for (unsigned int i = 1; i < numAccuracyBlock; i++) {
        m_accuracyMatrix[i - 1][m_curPhase - 1] /= (double) count;
    }
}



void MCMCAlg::initAccumuMargDistri() {
    unsigned i, j;
    for (i = 0; i < m_numVtx; i++) {
        for (j = 0; j < m_numTypeModel; j++) {
            m_accumuMargDistri[i][j] = 0.0;
        }
    }
}

void MCMCAlg::updateAccumuMargDistri(unsigned mutatevtxno, MCMC &mcmc) {
    const vector <pair<unsigned, double>> &lhvaripairs = mcmc.getLHVariPairs();
    for (unsigned int i = 0; i < lhvaripairs.size(); i++) {
        m_accumuMargDistri[mutatevtxno][lhvaripairs[i].first] += lhvaripairs[i].second;
    }
}

void MCMCAlg::checkTypeFixed() {
    if (!m_alltypefixed) {
        for (set<unsigned>::const_iterator csiiter = m_topVtxSet.begin(); csiiter != m_topVtxSet.end(); csiiter++) {
            unsigned int type = m_graph.getVertex(*csiiter).getType();
        }
    }
}

void MCMCAlg::updateBestTypeModelInPhase() {
    double bestLLHvalueInInit = max(m_MC_A->getBestLHvalue(), m_MC_B->getBestLHvalue());
    if (bestLLHvalueInInit > m_bestLLHvalueInPhase)
        m_bestLLHvalueInPhase = bestLLHvalueInInit;
    if (m_MC_A->getBestLHvalue() >= m_MC_B->getBestLHvalue()) {
        for (unsigned int i = 0; i < m_numVtx; i++)
            m_bestClfcInPhase[i] = m_MC_A->m_bestVtxTypeTable[i];
    } else {
        for (unsigned int i = 0; i < m_numVtx; i++)
            m_bestClfcInPhase[i] = m_MC_B->m_bestVtxTypeTable[i];
    }
}

