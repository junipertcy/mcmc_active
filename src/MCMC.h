#ifndef MCMC_H
#define MCMC_H

#include "TypeModel.h"

class MCMC {
public:
    explicit MCMC(TypeModel &tm, int blockmodeltype = 1, bool groupcorrected = false);

    unsigned getTargetType(unsigned int mutateVtxNo) noexcept;

    unsigned getMutateVtx() const { return m_mutateVtxNo; }

    void initVtxClassifiMatrix();

    void initGroupConnMatrix() { m_typeModel.initGroupConnMatrix(); }

    void updateVtxClassifiMatrix();

    void updateGroupConnMatrix() { m_typeModel.updateGroupConnMatrix(); }

    double **getVtxClassifiMatrix();

    void randInitTypeModel(const set<unsigned> &topVtxSet);

    void mutateTypeModel(unsigned v, unsigned t);

    void initBestTypeModel();

    void updateBestTypeModel();

    double getBestLHvalue() const { return m_bestLLHvalue; }

    const TypeModel &getTypeModel() const { return m_typeModel; }

    const vector<pair<unsigned, double>> &getLHVariPairs() const { return m_LHVariPairs; }

    double getLHvalue() { return m_logLikelihoodValue; }

    vector<pair<unsigned, double>> m_LHVariPairs;
    double_vec_t m_LLHVariTable;
    float_mat_t m_bestEdgeConnMatrix;
    float_vec_t m_bestGroupCardi;
    uint_vec_t m_bestVtxTypeTable;


private:
    void calcLHVari(unsigned int vtxNo, TypeModel &typeModel) noexcept;

    //case 2.1: undirected for model type 1
    void calcLHVariUDM1(unsigned int vtxNo, TypeModel &typeModel) noexcept;

    void initLogTable(unsigned int tablesize) noexcept;

    void initLogGammaTable(unsigned int tablesize) noexcept;

    unsigned calcTargetType() noexcept;

    double getLogFac(unsigned a);

    double getLog(unsigned a);

    double getLogDivFac(unsigned a, unsigned b);

    double calcLikelihood(const TypeModel &typemodel);

    double calcLikelihoodM1(const TypeModel &typemodel);

    //double calcGraphLikelihoodM2();
    void initLLHValue() {
        m_logLikelihoodValue = calcLikelihood(m_typeModel);
    }
    void updateLLHValue() { m_logLikelihoodValue += m_LLHVariTable[m_targetType]; }

    static double PI;
    double_vec_t m_transProbSelect;
    double_vec_t m_logfactable;
    double_vec_t m_logtable;
    double_vec_t m_loggammatable;  //log_gamma(0),log_gamma(0.5),log_gamma(1)...; log_gamma(i)=m_loggammatable[2i];

    unsigned int m_maxsizeLogtable;
    TypeModel &m_typeModel;

    unsigned int m_numVtx;
    unsigned int m_numType;

    unsigned int m_mutateVtxNo;
    unsigned int m_targetType;
    unsigned int m_sourceType;

    double MAXLOGDOUBLE;
    double MINLOGDOUBLE;
    static const unsigned OPTIMIZOR;

    double **dvtxClassifiMatrix;
    double m_bestLLHvalue;
    double m_logLikelihoodValue;
    unsigned int m_blockModelType;
    bool m_groupCorrected;

};

#endif


