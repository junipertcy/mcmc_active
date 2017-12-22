#ifndef MCMCALG_H
#define MCMCALG_H

#include "Learner.h"
#include "MCMC.h"

class MCMCAlg {
public:
    MCMCAlg(const Graph &graph, unsigned numTypeInModel, set<unsigned> &frozentypes, unsigned numOptInit,
            long numOptStep, unsigned numLearnerInit, long numLearnerStep, unsigned numPhase, unsigned numTop,
            unsigned learningMethod, unsigned modelType, bool groupcorrected);

    void runMCMCAlg();

    void runOnePhase();

    void runOneInit();

    void runOneStep();

private:
    int getRandRemainVtx();

    void getTopVtx();

    void initAccumuMargDistri();

    void initAccuracyMatrix();

    void updateAccuracyMatrix();

    void updateAccumuMargDistri(unsigned mutatevtxno, MCMC &mcmc);

    void checkTypeFixed();

    void updateBestTypeModelInPhase();

    static double stabRatio;
    static unsigned numAccuracyBlock;

    unsigned m_numVtx;
    unsigned m_numTypeModel;
    unsigned m_numTypeGraph;

    unsigned m_numOptInit;
    long m_numOptStep;
    unsigned m_numLearnerInit;
    long m_numLearnerStep;
    unsigned m_numPhase;
    unsigned m_numTop;
    unsigned m_method;

    set<unsigned> m_topVtxSet;
    set<unsigned> m_remainVtxSet;
    uint_vec_t m_topVtxSeq;

    const Graph &m_graph;

    std::unique_ptr<TypeModel> m_typeModelA;
    std::unique_ptr<TypeModel> m_typeModelB;
    std::unique_ptr<MCMC> m_MC_A;
    std::unique_ptr<MCMC> m_MC_B;
    std::unique_ptr<Learner> m_learner;

    long m_curStep;
    unsigned m_curPhase;
    unsigned m_curInit;

    unsigned m_resumePhase;

    float_mat_t m_accuracyMatrix;
    float_mat_t m_accumuMargDistri;
    uint_vec_t m_bestClfcInPhase;

    bool m_alltypefixed;

    double m_dSumLHvalueInPhase;

    long m_lNumLHvalueInPhase;

    double m_bestLLHvalueInPhase;


};

#endif
