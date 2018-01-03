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

    unsigned int get_Q() const noexcept;

    unsigned int getNumActiveType() const noexcept;

    unsigned int get_membership(unsigned int vtxno) const noexcept;

    bool isGTypeFrozen(unsigned int gtype) const noexcept;

    vector<set<unsigned int>> vtx_of_group_r_;
    uint_vec_t memberships_;
    
    uint_vec_t k_;
    uint_vec_t n_r_;
    uint_mat_t m_numTargetVtxGroup;
    uint_mat_t e_rs_;

    unsigned **m_numTargetGroupVtx;//m_numTargetGroupVtx[i][v] number of vertices in group i that has target of vertex v

private:
    unsigned int m_numActiveType;
    unsigned int N_;
    unsigned int Q_;
    set<unsigned int> m_frozenTypesInGraph;
    uu_map_t m_mapFzntyMty2Gty;
    uu_map_t m_mapFzntyGty2Mty;

    void randInitGroups(const set<unsigned int> &topVtxSet) noexcept;

    void mutate(unsigned int v, unsigned int t) noexcept;

};

#endif
