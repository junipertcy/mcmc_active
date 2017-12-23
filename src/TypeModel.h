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

    unsigned int getNumType() const noexcept { return m_numType; }

    unsigned int getNumActiveType() const noexcept { return m_numActiveType; }

    unsigned int getVtxType(unsigned int vtxno) const noexcept { return m_vtxTypeTable[vtxno]; }

    bool isGTypeFrozen(unsigned int gtype) const noexcept {
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
    unsigned int m_numActiveType;
    unsigned int m_numType;
    set<unsigned int> m_frozenTypesInGraph;
    uu_map_t m_mapFzntyMty2Gty;
    uu_map_t m_mapFzntyGty2Mty;

    void randInitGroups(const set<unsigned int> &topVtxSet) noexcept;

    void mutate(unsigned int v, unsigned int t) noexcept;

};

#endif
