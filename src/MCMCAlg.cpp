#include "MCMCAlg.h"
#include "Utility.cpp"

double MCMCAlg::stabRatio = 0.5;
unsigned int MCMCAlg::numAccuracyBlock = 10;

MCMCAlg::MCMCAlg(const Graph &graph, unsigned numTypeInModel, set<unsigned int> &frozentypes, unsigned int numOptInit,
                 long numOptStep, unsigned int numLearnerInit, long numLearnerStep, unsigned int numPhase, unsigned int numTop,
                 unsigned int learningMethod, unsigned int modelType, bool groupcorrected):
       m_graph(graph) {

    m_typeModelA = std::make_unique<TypeModel>(graph, numTypeInModel, frozentypes);
    m_typeModelB = std::make_unique<TypeModel>(graph, numTypeInModel, frozentypes);

    m_MC_A = std::make_unique<MCMC>(*m_typeModelA, modelType, groupcorrected);
    m_MC_B = std::make_unique<MCMC>(*m_typeModelB, modelType, groupcorrected);

    unsigned int i;
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

    N_ = m_typeModelA->getGraph().getNumVtx();
    m_numTypeGraph = m_typeModelA->getGraph().getNumType();
    m_numTypeModel = m_typeModelA->getNumType();
    for (i = 0; i < N_; i++) {
        m_remainVtxSet.insert(i);
    }
    //init top vtx seq
    for (i = 0; i < N_; i++) {
        unsigned gtype = graph.getVertex(i).getType();
        if (m_typeModelA->isGTypeFrozen(gtype)) {
            m_topVtxSet.insert(i);
            m_remainVtxSet.erase(i);
            m_topVtxSeq.push_back(i);
        }
    }
    //showiset(m_topVtxSet);
    m_bestClfcInPhase.resize(N_);
    m_accumuMargDistri.resize(N_);
    for (i = 0; i < N_; i++) {
        m_accumuMargDistri[i].resize(m_numTypeModel);
    }

    if (m_method == 3 || m_method == 4 || m_method == 5) {
        m_numLearnerInit = 0;
        m_numLearnerStep = 0;
    }

}

void MCMCAlg::runMCMCAlg() noexcept {
    /* IMPORTANT: Entry function */
    initAccuracyMatrix();
    for (m_curPhase_ = 1; (m_curPhase_ <= m_numPhase) && (!m_remainVtxSet.empty()); ++m_curPhase_) {
        runOnePhase();
    }
}

void MCMCAlg::initAccuracyMatrix() noexcept {
    m_accuracyMatrix.resize(numAccuracyBlock - 1);
    for (unsigned int i = 0; i < numAccuracyBlock - 1; i++)
        m_accuracyMatrix[i].resize(m_numPhase);
    for (unsigned int i = 0; i < numAccuracyBlock - 1; i++)
        for (unsigned int j = 0; j < m_numPhase; j++)
            m_accuracyMatrix[i][j] = 0.0;
}

void MCMCAlg::runOnePhase() noexcept {

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

    std::unique_ptr<double_mat_t> typeClassifiMatrix = std::make_unique<double_mat_t>();
    (*typeClassifiMatrix).resize(m_numTypeGraph);

    for (unsigned int i = 0; i < m_numTypeGraph; i++)
        (*typeClassifiMatrix)[i].resize(m_numTypeModel);
    for (unsigned int i = 0; i < m_numTypeGraph; i++) {
        for (unsigned int j = 0; j < m_numTypeModel; j++) {
            (*typeClassifiMatrix)[i][j] = 0.0;
        }
    }
    unsigned int realtype;
    double **vtxClassifiMatrixA = m_MC_A->getVtxClassifiMatrix();
    double **vtxClassifiMatrixB = m_MC_B->getVtxClassifiMatrix();
    for (unsigned int i = 0; i < N_; i++) {
        realtype = m_MC_A->getTypeModel().getGraph().getVertex(i).getType();
        for (unsigned int j = 0; j < m_numTypeModel; j++) {
            (*typeClassifiMatrix)[realtype][j] +=/*vtxClassifiMatrixA[i][j];//*/
                    (vtxClassifiMatrixA[i][j] + vtxClassifiMatrixB[i][j]) / 2;
        }
    }

    std::unique_ptr<uint_vec_t> topvtxGroupCardi = std::make_unique<uint_vec_t>();
    (*topvtxGroupCardi).resize(m_numTypeGraph);
    for (unsigned int i = 0; i < m_numTypeGraph; i++) {
        (*topvtxGroupCardi)[i] = 0;
    }

    for (auto const &i: m_topVtxSet) {
        realtype = (m_typeModelA->m_graph).getVertex(i).getType();
        (*topvtxGroupCardi)[realtype]++;
    }
}

void MCMCAlg::runOneInit() noexcept {


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

void MCMCAlg::runOneStep() noexcept {
    unsigned mutateVtxNo;
    unsigned sourceType;
    unsigned targetType;
    //typeModelA
    do {
        mutateVtxNo = (unsigned) getRandRemainVtx();
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
        mutateVtxNo = (unsigned) getRandRemainVtx();
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

    if ((m_curInit <= m_numLearnerInit) && (m_curStep > stabRatio * double(m_numLearnerStep)) &&
        (m_curStep <= m_numLearnerStep)) {
        m_learner->updateData();
    }
}

int MCMCAlg::getRandRemainVtx() noexcept {
    int index = randN((int)m_remainVtxSet.size());
    auto setiter = m_remainVtxSet.begin();
    for (unsigned int i = 0; i < index; ++i) {
        setiter++;
    }
    return *setiter;
}

void MCMCAlg::getTopVtx() noexcept {
    std::unique_ptr<uint_vec_t> arrayTop = std::make_unique<uint_vec_t>();
    (*arrayTop).resize(N_);

    unsigned int numtop = m_learner->getTopVtx(*arrayTop, m_numTop);
    for (unsigned int i = 0; i < numtop; i++) {
        m_topVtxSet.insert((*arrayTop)[i]);
        m_remainVtxSet.erase((*arrayTop)[i]);
        m_topVtxSeq.push_back((*arrayTop)[i]);
    }
}


void MCMCAlg::updateAccuracyMatrix() noexcept {
    unsigned int count = 0;
    double accuracyRatio = 1.0 / numAccuracyBlock;
    for (unsigned int vtxno = 0; vtxno < N_; vtxno++) {
        if (m_topVtxSet.count(vtxno)) {
            continue;
        }
        unsigned int realtype = m_typeModelA->getGraph().getVertex(vtxno).getType();
        double accuracyvalue = m_accumuMargDistri[vtxno][realtype];
        for (unsigned int i = 1; i < numAccuracyBlock; i++) {
            if (accuracyvalue >= accuracyRatio * i) {
                m_accuracyMatrix[i - 1][m_curPhase_ - 1]++;
            } else break;
        }
        count++;
    }
    for (unsigned int i = 1; i < numAccuracyBlock; i++) {
        m_accuracyMatrix[i - 1][m_curPhase_ - 1] /= (double) count;
    }
}



void MCMCAlg::initAccumuMargDistri() noexcept {
    for (unsigned int i = 0; i < N_; i++) {
        for (unsigned int j = 0; j < m_numTypeModel; j++) {
            m_accumuMargDistri[i][j] = 0.0;
        }
    }
}

void MCMCAlg::updateAccumuMargDistri(unsigned mutatevtxno, MCMC &mcmc) noexcept {
    const vector <pair<unsigned int, double>> &lhvaripairs = mcmc.getLHVariPairs();
    for (auto const &i: lhvaripairs) {
        m_accumuMargDistri[mutatevtxno][i.first] += i.second;
    }
}

void MCMCAlg::updateBestTypeModelInPhase() noexcept {
    double bestLLHvalueInInit = max(m_MC_A->getBestLHvalue(), m_MC_B->getBestLHvalue());
    if (bestLLHvalueInInit > m_bestLLHvalueInPhase)
        m_bestLLHvalueInPhase = bestLLHvalueInInit;
    if (m_MC_A->getBestLHvalue() >= m_MC_B->getBestLHvalue()) {
        for (unsigned int i = 0; i < N_; i++)
            m_bestClfcInPhase[i] = m_MC_A->m_bestVtxTypeTable[i];
    } else {
        for (unsigned int i = 0; i < N_; i++)
            m_bestClfcInPhase[i] = m_MC_B->m_bestVtxTypeTable[i];
    }
}

