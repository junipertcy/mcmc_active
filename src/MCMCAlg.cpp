#include "MCMC.h"
#include "Learner.h"

#include "MCMCAlg.h"
#include "Utility.cpp"

double MCMCAlg::stabRatio = 0.5;

MCMCAlg::MCMCAlg(const Graph &graph, unsigned numTypeInModel, set<unsigned int> &frozentypes, unsigned int numOptInit,
                 long numOptStep, unsigned int numLearnerInit, long numLearnerStep, unsigned int numPhase, unsigned int numTop,
                 unsigned int learningMethod, unsigned int modelType, bool groupcorrected):
       m_graph(graph)
{

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

    N_ = m_typeModelA->getGraph().get_N();
    Q_graph_ = m_typeModelA->getGraph().get_Q();
    Q_model_ = m_typeModelA->get_Q();
    for (i = 0; i < N_; i++) {
        m_remainVtxSet.insert(i);
    }

    //init top vtx seq
    for (i = 0; i < N_; i++) {
        unsigned int gtype = graph.getVertex(i).getType();
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
        m_learner->resetForNextPhase();
        for (m_curInit = 1; m_curInit <= max(m_numOptInit, m_numLearnerInit); m_curInit++) {
            runOneInit();
        }
        getTopVtx();
    }
}

void MCMCAlg::runOneInit() noexcept {
    m_MC_A->randInitTypeModel(this->m_topVtxSet);
    m_MC_B->randInitTypeModel(this->m_topVtxSet);

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
        mutateVtxNo = get_rand_remain_vtx();
        sourceType = m_typeModelA->memberships_[mutateVtxNo];
    } while (m_typeModelA->n_r_[sourceType] == 1);
    targetType = m_MC_A->getTargetType(mutateVtxNo);

    m_MC_A->gibbs_jump(mutateVtxNo, targetType);

    //typeModelB
    do {
        mutateVtxNo = get_rand_remain_vtx();
        sourceType = m_typeModelB->memberships_[mutateVtxNo];
    } while (m_typeModelB->n_r_[sourceType] == 1);
    targetType = m_MC_B->getTargetType(mutateVtxNo);

    m_MC_B->gibbs_jump(mutateVtxNo, targetType);

    if ((m_curInit <= m_numLearnerInit) && (m_curStep > stabRatio * double(m_numLearnerStep)) &&
        (m_curStep <= m_numLearnerStep)) {
        m_learner->updateData();
    }
}

unsigned int MCMCAlg::get_rand_remain_vtx() noexcept {
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

MCMCAlg::~MCMCAlg() {
    if (m_learner) {
        m_learner.release();
    }
}
