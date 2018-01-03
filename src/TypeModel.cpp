#include "TypeModel.h"
#include "Utility.cpp"
#include "output_functions.h"

TypeModel::TypeModel(const Graph &graph, unsigned numtype, set<unsigned> &frozentypes):
        m_graph(graph), Q_(numtype), N_(m_graph.get_N())
{
    unsigned int i;
    m_numActiveType = numtype - frozentypes.size();

    memberships_.resize(N_);
    n_r_.resize(Q_);
    k_.resize(Q_);

    m_numTargetVtxGroup.resize(N_);
    for (i = 0; i < N_; i++) {
        m_numTargetVtxGroup[i].resize(Q_);
    }

    m_numTargetGroupVtx = nullptr;
    e_rs_.resize(Q_);
    for (i = 0; i < Q_; i++) {
        e_rs_[i].resize(Q_);
    }

    auto suiciter = frozentypes.begin();
    for (i = m_numActiveType; suiciter != frozentypes.end(); suiciter++, i++) {
        m_mapFzntyMty2Gty.insert(uu_map_t::value_type(i, *suiciter));
        m_mapFzntyGty2Mty.insert(uu_map_t::value_type(*suiciter, i));
        m_frozenTypesInGraph.insert(*suiciter);
    }
}

unsigned int TypeModel::get_Q() const noexcept { return Q_; }

unsigned int TypeModel::getNumActiveType() const noexcept { return m_numActiveType; }

unsigned int TypeModel::get_membership(unsigned int vtxno) const noexcept { return memberships_[vtxno]; }

bool TypeModel::isGTypeFrozen(unsigned int gtype) const noexcept {
    if (m_frozenTypesInGraph.count(gtype) != 0) {
        return true;
    } else {
        return false;
    }
};

void TypeModel::randInitGroups(const set<unsigned> &topVtxSet) noexcept {
    unsigned int i, j;
    unsigned int type = 0;
    bool hasEmptyGroup;
    do {
        vtx_of_group_r_.clear();
        for (unsigned int q = 0; q < Q_; q++) {
            n_r_[q] = 0;
        }

        for (i = 0; i < N_; i++) {
            memberships_[i] = 0;
        }

        for (i = 0; i < Q_; i++) {
            set<unsigned int> groupSet;
            this->vtx_of_group_r_.push_back(groupSet);
        }

        for (i = 0; i < N_; i++) {
            if (topVtxSet.count(i)) {
                unsigned int gtype = m_graph.getVertex(i).getType();
                if (this->isGTypeFrozen(gtype)) {
                    type = m_mapFzntyGty2Mty[gtype];
                } else {
                    type = gtype;
                }
            } else {
                type = (unsigned) randN(m_numActiveType);
            }
            vtx_of_group_r_[type].insert(i);
            memberships_[i] = type;
            n_r_[type]++;
        }
        hasEmptyGroup = false;
        for (type = 0; type < m_numActiveType; type++) {
            if (this->n_r_[type] == 0) {
                hasEmptyGroup = true;
                break;
            }
        }
    } while (hasEmptyGroup);

    for (i = 0; i < N_; i++) {
        for (j = 0; j < Q_; j++) {
            m_numTargetVtxGroup[i][j] = 0;
        }
    }
    for (i = 0; i < N_; i++) {
        const set<unsigned int> &targets = m_graph.getVertex(i).getTargets();
        for (auto setiter: targets) {
            m_numTargetVtxGroup[i][memberships_[setiter]]++;
        }
    }

    if (m_numTargetGroupVtx) {
        for (i = 0; i < Q_; i++) {
            for (j = 0; j < N_; j++) {
                m_numTargetGroupVtx[i][j] = 0;
            }
        }
        for (i = 0; i < Q_; i++) {
            for (j = 0; j < N_; j++) {
                for (auto sourceVtx: vtx_of_group_r_[i]) {
                    if (m_graph.getVertex(sourceVtx).getTargets().count(j))
                        m_numTargetGroupVtx[i][j]++;
                }
            }
        }
    }

    for (i = 0; i < Q_; i++) {
        for (j = 0; j < Q_; j++) {
            e_rs_[i][j] = 0;
        }
    }

    for (i = 0; i < Q_; i++) {
        for (j = 0; j < Q_; j++) {
            for (auto sourceVtx: vtx_of_group_r_[i]) {
                e_rs_[i][j] += m_numTargetVtxGroup[sourceVtx][j];
            }
        }
    }
    for (i = 0; i < Q_; i++) {
        e_rs_[i][i] /= 2;
    }

    //group degrees data
    for (i = 0; i < Q_; i++) {
        k_[i] = 0;
        for (auto setuuiter: vtx_of_group_r_[i]) {
            k_[i] += m_graph.get_degree_at_v(setuuiter);
        }
    }
}

void TypeModel::mutate(unsigned int vtx, unsigned int target_group) noexcept {
    unsigned int source_group = memberships_[vtx];

    //update eij
    for (unsigned int q = 0; q < Q_; ++q) {
        if (q == source_group || q == target_group) {
            continue;
        }
        //eso'=eso-nvo, eos'=eso'
        e_rs_[source_group][q] -= m_numTargetVtxGroup[vtx][q];
        e_rs_[q][source_group] = e_rs_[source_group][q];
        //eto'=eto+nvo, eot'=eto'
        e_rs_[target_group][q] += m_numTargetVtxGroup[vtx][q];
        e_rs_[q][target_group] = e_rs_[target_group][q];
    }

    //est'=est-(nvt-nvs+lvv), ets'=est'
    e_rs_[source_group][target_group] -= (m_numTargetVtxGroup[vtx][target_group] - m_numTargetVtxGroup[vtx][source_group]);
    e_rs_[target_group][source_group] = e_rs_[source_group][target_group];

    //ess'=ess-nvs, ett'=ett+nvt+lvv;
    e_rs_[source_group][source_group] -= m_numTargetVtxGroup[vtx][source_group];
    e_rs_[target_group][target_group] += (m_numTargetVtxGroup[vtx][target_group]);

    //update nvi and niv
    const set<unsigned int> &vSources = m_graph.getVertex(vtx).getSources();
    for (auto const &i: vSources) {
        //nvs'=nvs-1, nvt'=nvt+1, here vtx is i;
        m_numTargetVtxGroup[i][source_group]--;
        m_numTargetVtxGroup[i][target_group]++;
    }

    //update ni
    n_r_[source_group]--;
    n_r_[target_group]++;
//    output_vec<uint_vec_t>(n_r_, std::clog);

    //update groupSets: source_group and target_group
    vtx_of_group_r_[source_group].erase(vtx);
    vtx_of_group_r_[target_group].insert(vtx);

    //update vtxTypeTable: vertex vtx has new type of target_group.
    memberships_[vtx] = target_group;

    //update group degrees
    k_[source_group] -= m_graph.get_degree_at_v(vtx);
    k_[target_group] += m_graph.get_degree_at_v(vtx);

}


