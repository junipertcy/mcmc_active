#ifndef MCMCALG_H
#define MCMCALG_H

#include "Learner.h"
#include "MCMC.h"

class MCMCAlg {
public:
    MCMCAlg(const Graph &graph, unsigned int numTypeInModel, set<unsigned int> &frozentypes, unsigned int numOptInit,
            long numOptStep, unsigned int numLearnerInit, long numLearnerStep, unsigned int numPhase, unsigned int numTop,
            unsigned int learningMethod, unsigned int modelType, bool groupcorrected);

    void runMCMCAlg() noexcept;

    void runOnePhase() noexcept;

    void runOneInit() noexcept;

    void runOneStep() noexcept;

private:
    unsigned int m_curPhase_;
    static unsigned int numAccuracyBlock;
    static double stabRatio;

    unsigned int N_;  // number of nodes of the network
    unsigned int Q_model_;
    unsigned int Q_graph_;

    int getRandRemainVtx() noexcept;

    void getTopVtx() noexcept;

    unsigned int m_numOptInit;
    long m_numOptStep;
    unsigned int m_numLearnerInit;
    long m_numLearnerStep;
    unsigned int m_numPhase;
    unsigned int m_numTop;
    unsigned int m_method;

    set<unsigned int> m_topVtxSet;
    set<unsigned int> m_remainVtxSet;

    uint_vec_t m_topVtxSeq;

    const Graph &m_graph;

    std::unique_ptr<TypeModel> m_typeModelA;
    std::unique_ptr<TypeModel> m_typeModelB;
    std::unique_ptr<MCMC> m_MC_A;
    std::unique_ptr<MCMC> m_MC_B;
    std::unique_ptr<Learner> m_learner;

    long m_curStep;

    unsigned int m_curInit;
};

#endif
