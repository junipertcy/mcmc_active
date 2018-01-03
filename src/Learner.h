#ifndef LEARNER_H
#define LEARNER_H

#include "MCMC.h"

bool compare(const pair<double,unsigned int>&i, const pair<double,unsigned int>&j) noexcept;

class Learner {

public:
    virtual void resetForNextPhase() noexcept { ; }

    virtual void resetForNextInit() noexcept { ; }

    virtual unsigned int getTopVtx(uint_vec_t &arrayTop, unsigned int numTop) noexcept { return 0; }

    virtual void updateData() noexcept { ; }

    virtual unsigned int get_N() const noexcept { return m_MC_A.getTypeModel().getGraph().get_N(); }

    virtual unsigned int getNumTypeGraph() const noexcept { return m_MC_A.getTypeModel().getGraph().get_Q(); }

    virtual unsigned int getNumTypeModel() const noexcept { return m_MC_A.getTypeModel().get_Q(); }

protected:
    Learner(const MCMC &mca, const MCMC &mcb, const set<unsigned> &topVtxSet);

    const MCMC &m_MC_A;
    const MCMC &m_MC_B;
//// TODO: it was originally declared as const; how to retain the type status?
//    std::shared_ptr<MCMC> m_MC_A;
//    std::shared_ptr<MCMC> m_MC_B;

    const set<unsigned int> &m_topVtxSet;
    unsigned int N_;// number of vertices in graph
    unsigned int Q_;// number of types in the model (may be not that in the real graph)

};

class MutualInfo : public Learner {
public:
    MutualInfo(const MCMC &mca, const MCMC &mcb, const set<unsigned> &topVtxSet);

    void resetForNextPhase() noexcept override;

    void resetForNextInit() noexcept override;

    unsigned int getTopVtx(uint_vec_t &arrayTop, unsigned numTop) noexcept override;

    void updateData() noexcept override;

    unsigned int get_N() const noexcept override;

    unsigned int getNumTypeGraph() const noexcept override;

    unsigned int getNumTypeModel() const noexcept override;

private:
    double_mat_t m_accumuMargDistri;
    float_vec_t m_accumuCondEntropy;
    float_vec_t m_numAccumuCondEntropy;
};


#endif