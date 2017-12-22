#ifndef LEARNER_H
#define LEARNER_H

#include "MCMC.h"

class Learner {

public:
    virtual void resetForNextPhase() noexcept { return; }

    virtual void resetForNextInit() noexcept { return; }

    virtual unsigned int getTopVtx(unsigned int *arrayTop, unsigned numTop) noexcept { return 0; }

    virtual void updateData() noexcept { return; }

    virtual unsigned int getNumVtx() const noexcept { return m_MC_A.getTypeModel().getGraph().getNumVtx(); }

    virtual unsigned int getNumTypeGraph() const noexcept { return m_MC_A.getTypeModel().getGraph().getNumType(); }

    virtual unsigned int getNumTypeModel() const noexcept { return m_MC_A.getTypeModel().getNumType(); }

protected:
    Learner(const MCMC &mca, const MCMC &mcb, const set<unsigned> &topVtxSet, string methodname);

    const MCMC &m_MC_A;
    const MCMC &m_MC_B;
    const set<unsigned> &m_topVtxSet;
    unsigned m_numVtx;// number of vertices in graph
    unsigned m_numType;// number of types in the model (may be not that in the real graph)
    static const unsigned QUERYSTRATEGY;// query strategy
    string m_methodName;

    unsigned m_arraysize;
    unsigned m_numtop;

    uint_vec_t m_arrayTopVtxNo;  //some query strategies may output more than "numTop" vertices
    uint_vec_t m_arrayVtxNo;
    float_vec_t m_arrayLearnScores;

    unsigned int *m_arrayVtxNoSort;
    float *m_arrayLearnScoresSort;

};

class MutualInfo : public Learner {
public:
    MutualInfo(const MCMC &mca, const MCMC &mcb, const set<unsigned> &topVtxSet);

    void resetForNextPhase() noexcept override;

    void resetForNextInit() noexcept override;

    unsigned int getTopVtx(unsigned int *arrayTop, unsigned numTop) noexcept override;

    void updateData() noexcept override;

    unsigned int getNumVtx() const noexcept override;

    unsigned int getNumTypeGraph() const noexcept override;

    unsigned int getNumTypeModel() const noexcept override;

private:
    float_mat_t m_accumuMargDistri;
    float_vec_t m_accumuCondEntropy;
    float_vec_t m_numAccumuCondEntropy;

    unsigned int updateDataIndex;
};






#endif