#include "Graph.h"


Graph::Graph(string sGraphFileName) {
    unsigned i;
    ifstream infile(sGraphFileName.c_str());
    string words;

    if (!infile) {
        cerr << "unable to open the input graph data file:"
             << sGraphFileName << endl;
    }
    set<unsigned> selfloopVtxSet;
    m_numVtx = 0;
    m_numType = 0;
    set <string> seenIds;
    set <string> seenValues;
    string word;
    string id;
    string value;
    unsigned type;
    string sourceid;
    string targetid;
    unsigned sourceindex;
    unsigned targetindex;
    setSelfloop(false);
    setDirected(false);
    while (infile >> word) {
        if (word == "directed") {
            infile >> word;
            if (word == "1") {
                this->setDirected(true);
            } else {
                this->setDirected(false);
            }
            break;
        }
    }
    while (infile >> word) {
        if (word == "node") {

            Vertex vtx = Vertex();
            while (infile >> word) {
                if (word == "id") {
                    infile >> id;
                    //id=atoi(word.c_str());
                    vtx.setId(id);
                }
                if (word == "type" || word == "value") {
                    infile >> value;
                    vtx.setValue(value);
                }
                if (word == "]")
                    break;
            }
            /************************************/
            if (seenIds.count(id)) {
                continue;
            } else {
                vtx.setIndex(m_numVtx);
                m_mapId2Index.insert(valType_su(id, m_numVtx));
                seenIds.insert(id);
            }
            if (seenValues.count(value)) {
                type = m_mapValue2Type[value];
                vtx.setType(type);
            } else {
                type = m_numType;
                vtx.setType(m_numType);
                m_mapValue2Type.insert(valType_su(value, m_numType));
                seenValues.insert(value);
                m_numType++;
            }
            m_mapIndex2Id.insert(valType_us(m_numVtx, id));
            m_mapType2Value.insert(valType_us(type, value));
            m_numVtx++;
            this->vtxList.push_back(vtx);
            /*************************/

        }
        if (word == "edge") {
            while (infile >> word) {
                if (word == "source") {
                    infile >> sourceid;
                }
                if (word == "target") {
                    infile >> targetid;
                }
                if (word == "]")
                    break;
            }
            if (m_mapId2Index.count(sourceid)) {
                sourceindex = m_mapId2Index[sourceid];
            } else {
                cerr << "unable to find the id information of the source vertex" << endl;
                continue;
            }
            if (m_mapId2Index.count(targetid)) {
                targetindex = m_mapId2Index[targetid];
            } else {
                cerr << "unable to find the id information of the target vertex" << endl;
                continue;
            }
            if (sourceid == targetid) {
                selfloopVtxSet.insert(sourceindex);
                setSelfloop(true);
            }
            vtxList[sourceindex].addTarget(targetindex);
            vtxList[targetindex].addSource(sourceindex);
            vtxList[sourceindex].addSource(targetindex);
            vtxList[targetindex].addTarget(sourceindex);
        }
    }

    infile.close();

    vtxSelfloopFlag.assign(m_numVtx, 0);
    if (m_selfloop) {
        for (i = 0; i < m_numVtx; i++) {
            if (selfloopVtxSet.count(i))
                vtxSelfloopFlag[i] = 1;
        }
    }

    m_groupConnNumMatrix.resize(m_numType);
    for (i = 0; i < m_numType; i++) {
        m_groupConnNumMatrix[i].resize(m_numType);
    }
    m_groupCardi.resize(m_numType);
    m_groupConnMatrix.resize(m_numType);

    for (i = 0; i < m_numType; i++) {
        m_groupConnMatrix[i].resize(m_numType);
    }

    calcNumEdges();

    calcGroupConnNumMatrix();
    calcGroupConnMatrix();
}

void Graph::calcNumEdges() {
    unsigned int numEgs = 0;
    for (auto const &i: vtxList) {
        numEgs += i.getOutDegree();
    }
    m_numEgs = numEgs / 2;
}

void Graph::calcGroupConnNumMatrix() {
    unsigned int i, j;
    unsigned int stype, ttype;

    for (i = 0; i < m_numType; i++) {
        for (j = 0; j < m_numType; j++) {
            m_groupConnNumMatrix[i][j] = 0;
        }
    }

    for (i = 0; i < m_numType; i++) {
        m_groupCardi[i] = 0;
    }

    for (i = 0; i < vtxList.size(); i++) {
        const set<unsigned> &edges = vtxList[i].getTargets();

        stype = vtxList[i].getType();
        m_groupCardi[stype]++;

        for (auto const &k: edges) {
            ttype = vtxList[k].getType();
            if ((!isDirected()) && (ttype == stype) && (k > i))
                continue;
            m_groupConnNumMatrix[stype][ttype]++;
        }
    }
}

void Graph::calcGroupConnMatrix() {
    unsigned i, j;
    for (i = 0; i < m_numType; i++) {
        for (j = i; j < m_numType; j++) {
            if (m_groupConnNumMatrix[i][j] == 0) {
                m_groupConnMatrix[i][j] = 0.0;
                continue;
            }
            if (i == j) {
                m_groupConnMatrix[i][i] = ((double) (2 * m_groupConnNumMatrix[i][i])) /
                                          (((double) m_groupCardi[i]) *
                                           ((double) (m_groupCardi[i] - 1))); //  2eii/(ni*(ni-1))
            } else {
                m_groupConnMatrix[i][j] = ((double) m_groupConnNumMatrix[i][j]) / (((double) m_groupCardi[i]) *
                                                                                   ((double) (m_groupCardi[j]))); //  eij/(ni*nj)
            }
        }
    }
    for (i = 1; i < m_numType; i++) {
        for (j = 0; j < i; j++) {
            m_groupConnMatrix[i][j] = m_groupConnMatrix[j][i];
        }
    }
}

const Graph::Vertex &Graph::getVertex(unsigned vtxno) const {
    if (vtxno < m_numVtx)
        return vtxList[vtxno];
    else {
        cerr << "Cannot get the vertex, the vertex no is too large." << endl
             << "-- getVertex()::Graph" << endl;
        return vtxList[m_numVtx - 1];
    }
}

const unsigned Graph::getVtxDegree(unsigned vtxno) const {
    if (vtxno > vtxList.size()) {
        cerr << "the vtx no provided is illegal.-- Graph::getVtxDegree()" << endl;
        return 0;
    }
    if (vtxSelfloopFlag[vtxno] == 0) {
        return vtxList[vtxno].getInDegree();
    } else {
        return vtxList[vtxno].getInDegree() + 1;
    }
}