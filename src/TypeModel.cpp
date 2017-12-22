#include "TypeModel.h"
#include "Utility.cpp"

TypeModel::TypeModel(const Graph &graph, unsigned numtype, set<unsigned> &frozentypes):
        m_graph(graph), m_numType(numtype)
{
    unsigned i;
    m_numActiveType = numtype - frozentypes.size();
    unsigned int m_numVtx = m_graph.getNumVtx();

    m_vtxTypeTable.resize(m_numVtx);
    m_groupCardiTable.resize(m_numVtx);
    m_groupDegrees.resize(m_numType);

    m_numTargetVtxGroup.resize(m_numVtx);
    for (i = 0; i < m_numVtx; i++) {
        m_numTargetVtxGroup[i].resize(m_numType);
    }

    m_numTargetGroupVtx = 0;
    m_numEdgesOf2Groups.resize(m_numType);
    for (i = 0; i < m_numType; i++) {
        m_numEdgesOf2Groups[i].resize(m_numType);
    }

    groupConnMatrix.resize(m_numType);
    for (i = 0; i < m_numType; i++) {
        groupConnMatrix[i].resize(m_numType);
    }

    lvtxClassifiMatrix.resize(m_numVtx);
    dvtxClassifiMatrix.resize(m_numVtx);
    for (i = 0; i < m_numVtx; i++) {
        lvtxClassifiMatrix[i].resize(m_numType);
        dvtxClassifiMatrix[i].resize(m_numType);
    }

    auto suiciter = frozentypes.begin();
    for (i = m_numActiveType; suiciter != frozentypes.end(); suiciter++, i++) {
        m_mapFzntyMty2Gty.insert(valType_uu(i, *suiciter));
        m_mapFzntyGty2Mty.insert(valType_uu(*suiciter, i));
        m_frozenTypesInGraph.insert(*suiciter);
    }
}

void TypeModel::randInitGroups(const set<unsigned> &topVtxSet) {
    unsigned i, j;
    unsigned type = 0;
    unsigned int m_numVtx = m_graph.getNumVtx();
    bool hasEmptyGroup;
    do {
        m_groupSets.clear();
        for (i = 0; i < m_numType; i++) {
            m_groupCardiTable[i] = 0;
        }

        for (i = 0; i < m_numVtx; i++) {
            m_vtxTypeTable[i] = 0;
        }

        for (i = 0; i < m_numType; i++) {
            set<unsigned> groupSet;
            this->m_groupSets.push_back(groupSet);
        }

        for (i = 0; i < m_numVtx; i++) {
            if (topVtxSet.count(i)) {
                unsigned gtype = m_graph.getVertex(i).getType();
                if (this->isGTypeFrozen(gtype))
                    type = m_mapFzntyGty2Mty[gtype];
                else
                    type = gtype;
            } else {
                type = (unsigned) randN(m_numActiveType);
            }
            m_groupSets[type].insert(i);
            m_vtxTypeTable[i] = type;
            m_groupCardiTable[type]++;
        }
        hasEmptyGroup = false;
        for (type = 0; type < m_numActiveType; type++) {
            if (this->m_groupCardiTable[type] == 0) {
                hasEmptyGroup = true;
                break;
            }
        }
    } while (hasEmptyGroup);

    for (i = 0; i < m_numVtx; i++) {
        for (j = 0; j < m_numType; j++) {
            m_numTargetVtxGroup[i][j] = 0;
        }
    }
    set<unsigned>::const_iterator setiter;
    for (i = 0; i < m_numVtx; i++) {
        const set<unsigned> &targets = m_graph.getVertex(i).getTargets();
        for (setiter = targets.begin(); setiter != targets.end(); setiter++)
            m_numTargetVtxGroup[i][m_vtxTypeTable[*setiter]]++;
    }

    set<unsigned>::const_iterator siiter;
    if (m_numTargetGroupVtx) {
        for (i = 0; i < m_numType; i++) {
            for (j = 0; j < m_numVtx; j++) {
                m_numTargetGroupVtx[i][j] = 0;
            }
        }
        for (i = 0; i < m_numType; i++) {
            const set<unsigned> &groupSet = m_groupSets[i];
            for (j = 0; j < m_numVtx; j++) {
                for (siiter = groupSet.begin(); siiter != groupSet.end(); siiter++) {
                    unsigned sourceVtx = *siiter;
                    if (m_graph.getVertex(sourceVtx).getTargets().count(j))
                        m_numTargetGroupVtx[i][j]++;
                }
            }
        }
    }
    for (i = 0; i < m_numType; i++) {
        for (j = 0; j < m_numType; j++) {
            m_numEdgesOf2Groups[i][j] = 0;
        }
    }
    for (i = 0; i < m_numType; i++) {
        for (j = 0; j < m_numType; j++) {
            const set<unsigned> &groupSet = m_groupSets[i];
            for (siiter = groupSet.begin(); siiter != groupSet.end(); siiter++) {
                unsigned sourceVtx = *siiter;
                m_numEdgesOf2Groups[i][j] += m_numTargetVtxGroup[sourceVtx][j];
            }
        }
    }
    for (i = 0; i < m_numType; i++) {
        unsigned numSelfloop = 0;
        const set<unsigned> &groupSet = m_groupSets[i];
        for (siiter = groupSet.begin(); siiter != groupSet.end(); siiter++) {
            unsigned sourceVtx = *siiter;
            if (m_graph.vtxHasSelfloop(sourceVtx))
                numSelfloop++;
        }
        m_numEdgesOf2Groups[i][i] = (m_numEdgesOf2Groups[i][i] + numSelfloop) / 2;
    }

    //group degrees data
    set<unsigned>::const_iterator setuuiter;
    for (i = 0; i < m_numType; i++) {
        setuuiter = m_groupSets[i].begin();
        m_groupDegrees[i] = 0;
        for (; setuuiter != m_groupSets[i].end(); setuuiter++) {
            m_groupDegrees[i] += m_graph.getVtxDegree(*setuuiter);
        }
    }
}

void TypeModel::mutate(unsigned v, unsigned t) {
    unsigned o;
    unsigned s = m_vtxTypeTable[v];

    unsigned lvv = 0;
    if (m_graph.vtxHasSelfloop(v))
        lvv = 1;
    //update eij
    for (o = 0; o < m_numType; o++) {
        if (o == s || o == t)
            continue;
        //eso'=eso-nvo, eos'=eso'
        m_numEdgesOf2Groups[s][o] -= m_numTargetVtxGroup[v][o];
        m_numEdgesOf2Groups[o][s] = m_numEdgesOf2Groups[s][o];
        //eto'=eto+nvo, eot'=eto'
        m_numEdgesOf2Groups[t][o] += m_numTargetVtxGroup[v][o];
        m_numEdgesOf2Groups[o][t] = m_numEdgesOf2Groups[t][o];
    }
    //est'=est-(nvt-nvs+lvv), ets'=est'
    m_numEdgesOf2Groups[s][t] -= (m_numTargetVtxGroup[v][t] - m_numTargetVtxGroup[v][s] + lvv);
    m_numEdgesOf2Groups[t][s] = m_numEdgesOf2Groups[s][t];
    //ess'=ess-nvs, ett'=ett+nvt+lvv;
    m_numEdgesOf2Groups[s][s] -= m_numTargetVtxGroup[v][s];
    m_numEdgesOf2Groups[t][t] += (m_numTargetVtxGroup[v][t] + lvv);

    //update nvi and niv
    const set<unsigned> &vSources = m_graph.getVertex(v).getSources();
    set<unsigned>::const_iterator setiter = vSources.begin();
    for (; setiter != vSources.end(); setiter++) {
        //nvs'=nvs-1, nvt'=nvt+1, here v is i;
        m_numTargetVtxGroup[*setiter][s]--;
        m_numTargetVtxGroup[*setiter][t]++;
    }

    //update ni
    m_groupCardiTable[s]--;
    m_groupCardiTable[t]++;
    //update groupSets: s and t
    m_groupSets[s].erase(v);
    m_groupSets[t].insert(v);
    //update vtxTypeTable: vertex v has new type of t.
    m_vtxTypeTable[v] = t;
    //update group degrees
    m_groupDegrees[s] -= m_graph.getVtxDegree(v);
    m_groupDegrees[t] += m_graph.getVtxDegree(v);

}

void TypeModel::initGroupConnMatrix() {
    unsigned i, j;
    for (i = 0; i < m_numType; i++) {
        for (j = 0; j < m_numType; j++) {
            groupConnMatrix[i][j] = 0;
        }
    }
    numAccuGCM = 0l;
}

void TypeModel::updateGroupConnMatrix() {
    unsigned i, j;
    for (i = 0; i < m_numType; i++) {
        for (j = i; j < m_numType; j++) {
            if (m_numEdgesOf2Groups[i][j] == 0) {
                continue;
            }
            if (i == j)
                groupConnMatrix[i][i] += ((double) (2 * m_numEdgesOf2Groups[i][i])) /
                                         (((double) m_groupCardiTable[i]) *
                                          ((double) (m_groupCardiTable[i] - 1))); //  2eii/(ni*(ni-1))
            else
                groupConnMatrix[i][j] += ((double) m_numEdgesOf2Groups[i][j]) / (((double) m_groupCardiTable[i]) *
                                                                                 ((double) (m_groupCardiTable[j]))); //  eij/(ni*nj)
        }
    }
    for (i = 1; i < m_numType; i++) {
        for (j = 0; j < i; j++) {
            groupConnMatrix[i][j] = groupConnMatrix[j][i];
        }
    }
    numAccuGCM++;
}

