#ifndef TYPEMODEL_H
#define TYPEMODEL_H

#include "Graph.h"
#include <vector>
#include <set>

class TypeModel {
    friend class MCMC;

public:
    const Graph &m_graph;

    const Graph &getGraph() const noexcept { return m_graph; }

    TypeModel(const Graph &graph, unsigned numtype/*including frozen types*/,
              set<unsigned> &frozentypes/*frozen types defined in graph*/);

    unsigned getNumType() const { return m_numType; }

    unsigned getNumActiveType() const { return m_numActiveType; }

    unsigned getVtxType(unsigned vtxno) const { return m_vtxTypeTable[vtxno]; }

    bool isGTypeFrozen(unsigned gtype) const {
        if (m_frozenTypesInGraph.count(gtype) != 0)
            return true;
        else
            return false;
    };

    vector<set<unsigned int>> m_groupSets;
    uint_vec_t m_vtxTypeTable;
    uint_vec_t m_groupDegrees;
    uint_vec_t m_groupCardiTable;
    uint_mat_t m_numTargetVtxGroup;
    uint_mat_t m_numEdgesOf2Groups;

    unsigned **m_numTargetGroupVtx;//m_numTargetGroupVtx[i][v] number of vertices in group i that has target of vertex v

private:
    uint_mat_t groupConnMatrix;
    float_mat_t lvtxClassifiMatrix;
    float_mat_t dvtxClassifiMatrix;

    void initGroupConnMatrix();

    void updateGroupConnMatrix();

    void randInitGroups(const set<unsigned int> &topVtxSet);

    void mutate(unsigned int v, unsigned int t);

    unsigned int m_numActiveType;
    unsigned int m_numType;
    set<unsigned> m_frozenTypesInGraph;
    map<unsigned, unsigned> m_mapFzntyMty2Gty;
    map<unsigned, unsigned> m_mapFzntyGty2Mty;
    long long numAccuGCM;
};

#endif
