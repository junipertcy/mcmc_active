#include "MCMCAlg.h"
#include "Utility.cpp"
#include "output_functions.h"

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
    Q_graph_ = m_typeModelA->getGraph().getNumType();
    Q_model_ = m_typeModelA->getNumType();
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

}

void MCMCAlg::runMCMCAlg() noexcept {
    /* IMPORTANT: Entry function */
    for (m_curPhase_ = 1; (m_curPhase_ <= m_numPhase) && (!m_remainVtxSet.empty()); ++m_curPhase_) {
        runOnePhase();
    }
}


void MCMCAlg::runOnePhase() noexcept {
    m_learner->resetForNextPhase();

    (this->m_MC_A)->initVtxClassifiMatrix();
    (this->m_MC_B)->initVtxClassifiMatrix();

    for (m_curInit = 1; m_curInit <= max(m_numOptInit, m_numLearnerInit); m_curInit++) {
        runOneInit();
    }

    getTopVtx();

    std::unique_ptr<double_mat_t> typeClassifiMatrix = std::make_unique<double_mat_t>();
    (*typeClassifiMatrix).resize(Q_graph_);

    for (unsigned int i = 0; i < Q_graph_; i++)
        (*typeClassifiMatrix)[i].resize(Q_model_);
    for (unsigned int i = 0; i < Q_graph_; i++) {
        for (unsigned int j = 0; j < Q_model_; j++) {
            (*typeClassifiMatrix)[i][j] = 0.0;
        }
    }
    unsigned int realtype;
    double **vtxClassifiMatrixA = m_MC_A->getVtxClassifiMatrix();
    double **vtxClassifiMatrixB = m_MC_B->getVtxClassifiMatrix();
    for (unsigned int i = 0; i < N_; i++) {
        realtype = m_MC_A->getTypeModel().getGraph().getVertex(i).getType();
        for (unsigned int j = 0; j < Q_model_; j++) {
            (*typeClassifiMatrix)[realtype][j] +=/*vtxClassifiMatrixA[i][j];//*/
                    (vtxClassifiMatrixA[i][j] + vtxClassifiMatrixB[i][j]) / 2;
        }
    }

    std::unique_ptr<uint_vec_t> topvtxGroupCardi = std::make_unique<uint_vec_t>();
    (*topvtxGroupCardi).resize(Q_graph_);
    for (unsigned int i = 0; i < Q_graph_; i++) {
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
        if (m_numOptInit < m_numLearnerInit) {
            maxstep = m_numLearnerStep;
        } else {
            maxstep = m_numOptStep;
        }
    }

    for (m_curStep = 1; m_curStep <= maxstep; m_curStep++) {
        runOneStep();
    }
}

void MCMCAlg::runOneStep() noexcept {
    unsigned int mutateVtxNo;
    unsigned int sourceType;
    unsigned int targetType;
    //typeModelA
    do {
        mutateVtxNo = (unsigned) getRandRemainVtx();
        sourceType = m_typeModelA->memberships_[mutateVtxNo];
    } while (m_typeModelA->m_groupCardiTable[sourceType] == 1);
    targetType = m_MC_A->getTargetType(mutateVtxNo);

    m_MC_A->mutateTypeModel(mutateVtxNo, targetType);


    if ((m_curInit <= m_numOptInit) && (m_curStep > stabRatio * double(m_numOptStep)) && (m_curStep <= m_numOptStep)) {
        (this->m_MC_A)->updateVtxClassifiMatrix();
        (this->m_MC_A)->updateBestTypeModel();
    }

    //typeModelB
    do {
        mutateVtxNo = (unsigned) getRandRemainVtx();
        sourceType = m_typeModelB->memberships_[mutateVtxNo];
    } while (m_typeModelB->m_groupCardiTable[sourceType] == 1);
    targetType = m_MC_B->getTargetType(mutateVtxNo);

    m_MC_B->mutateTypeModel(mutateVtxNo, targetType);

    if ((m_curInit <= m_numOptInit) && (m_curStep > stabRatio * double(m_numOptStep)) && (m_curStep <= m_numOptStep)) {
        (this->m_MC_B)->updateVtxClassifiMatrix();
        (this->m_MC_B)->updateBestTypeModel();
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