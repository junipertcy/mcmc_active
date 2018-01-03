#ifndef MCMC_H
#define MCMC_H

#include "TypeModel.h"

class MCMC {
public:
    explicit MCMC(TypeModel &tm, int blockmodeltype = 1, bool groupcorrected = false);

    unsigned int getTargetType(unsigned int mutateVtxNo) noexcept;

    unsigned int get_selected_vtx() const noexcept;

    void randInitTypeModel(const set<unsigned int> &topVtxSet) noexcept;

    void gibbs_jump(unsigned int v, unsigned int t) noexcept;

    const TypeModel &getTypeModel() const noexcept;

    const vector<pair<unsigned int, double>> &get_likelihood_variation_pairs() const noexcept;

    vector<pair<unsigned int, double>> m_LHVariPairs;
    double_vec_t m_LLHVariTable;


private:
    void compute_likelihood(unsigned int vtxNo, TypeModel &typeModel) noexcept;

    void _compute_likelihood(unsigned int vtxNo, TypeModel &typeModel) noexcept;

    void initLogTable(unsigned int tablesize) noexcept;

    void initLogGammaTable(unsigned int tablesize) noexcept;

    unsigned int compute_gibbs_jump() noexcept;

    double getLogFac(unsigned int a) noexcept;

    double getLog(unsigned int a) noexcept;

    double getLogDivFac(unsigned int a, unsigned int b) noexcept;

    double calcLikelihood(const TypeModel &typemodel) noexcept;

    double calcLikelihoodM1(const TypeModel &typemodel) noexcept;

    //double calcGraphLikelihoodM2();
    void init_log_likelihood() noexcept;

    void update_log_likelihood() noexcept;

    static double PI;
    double_vec_t m_transProbSelect;
    double_vec_t m_logfactable;
    double_vec_t m_logtable;
    double_vec_t m_loggammatable;  //log_gamma(0),log_gamma(0.5),log_gamma(1)...; log_gamma(i)=m_loggammatable[2i];

    unsigned int m_maxsizeLogtable;
    TypeModel &m_typeModel;

    unsigned int N_;
    unsigned int Q_;

    unsigned int m_mutateVtxNo;
    unsigned int m_targetType;
    unsigned int m_sourceType;

    double MAXLOGDOUBLE;
    double MINLOGDOUBLE;
    static const unsigned OPTIMIZOR;

    double **dvtxClassifiMatrix;
    double m_logLikelihoodValue;
    unsigned int m_blockModelType;
    bool m_groupCorrected;
};

#endif


