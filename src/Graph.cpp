#include "Graph.h"


Graph::Graph(string sGraphFileName) {
    unsigned i;
    bool isblogud = false;
    ifstream infile(sGraphFileName.c_str());

    /************************************/
    string vtxdeg0filename("vtxDeg0IDs.txt");
    ifstream vtxdeg0file(vtxdeg0filename.c_str());
    set <string> vtxdeg0ids;
    string words;
    /************************************/

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
        if (word.compare("directed") == 0) {
            infile >> word;
            if (word.compare("1") == 0)
                this->setDirected(true);
            else
                this->setDirected(false);
            break;
        }
    }
    while (infile >> word) {
        if (word.compare("node") == 0) {

            Vertex vtx = Vertex();
            while (infile >> word) {
                if (word.compare("id") == 0) {
                    infile >> id;
                    //id=atoi(word.c_str());
                    vtx.setId(id);
                }
                if (word.compare("type") == 0 || word.compare("value") == 0) {
                    infile >> value;
                    vtx.setValue(value);
                }
                if (word.compare("]") == 0)
                    break;
            }
            /************************************/
            if (isblogud && (vtxdeg0ids.count(id))) { ;
            } else {
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
            }
            /*************************/

        }
        if (word.compare("edge") == 0) {
            while (infile >> word) {
                if (word.compare("source") == 0) {
                    infile >> sourceid;
                }
                if (word.compare("target") == 0) {
                    infile >> targetid;
                }
                if (word.compare("]") == 0)
                    break;
            }
            if (m_mapId2Index.count(sourceid)) {
                sourceindex = m_mapId2Index[sourceid];
            } else {
                sourceindex = 0;
                cerr << "unable to find the id information of the source vertex" << endl;
                continue;
            }
            if (m_mapId2Index.count(targetid)) {
                targetindex = m_mapId2Index[targetid];
            } else {
                targetindex = 0;
                cerr << "unable to find the id information of the target vertex" << endl;
                continue;
            }
            if (sourceid.compare(targetid) == 0) {
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

    m_numSelfloop = selfloopVtxSet.size();

    vtxSelfloopFlag.assign(m_numVtx, 0);
    if (m_selfloop) {
        for (i = 0; i < m_numVtx; i++) {
            if (selfloopVtxSet.count(i))
                vtxSelfloopFlag[i] = 1;
        }
    }

    //
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
//
//Graph::~Graph(void) {
//    unsigned i;
//    for (i = 0; i < m_numType; i++) {
//        delete[] m_groupConnNumMatrix[i];
//    }
//    delete[] m_groupConnNumMatrix;
//
//    for (i = 0; i < m_numType; i++) {
//        delete[] m_groupConnMatrix[i];
//    }
//    delete[] m_groupConnMatrix;
//    delete[] m_groupCardi;
//}


void Graph::calcNumEdges() {
    unsigned int numEgs = 0;
    for (unsigned int i = 0; i < vtxList.size(); ++i) {
        numEgs += vtxList[i].getOutDegree();
    }
    m_numEgs = numEgs / 2;
}

void Graph::calcGroupConnNumMatrix() {
    unsigned i, j;
    unsigned stype, ttype;

    for (i = 0; i < m_numType; i++) {
        for (j = 0; j < m_numType; j++) {
            m_groupConnNumMatrix[i][j] = 0;
        }
    }
    for (i = 0; i < m_numType; i++)
        m_groupCardi[i] = 0;
    //unsigned numNodesD0=0;
    for (i = 0; i < vtxList.size(); i++) {
        const set<unsigned> &edges = vtxList[i].getTargets();
        set<unsigned>::const_iterator siter = edges.begin();
        stype = vtxList[i].getType();
        m_groupCardi[stype]++;
        for (; siter != edges.end(); siter++) {
            ttype = vtxList[*siter].getType();
            if ((!isDirected()) && (ttype == stype) && ((*siter) > i))
                continue;
            m_groupConnNumMatrix[stype][ttype]++;
        }
    }
    //cout<<"Number of Nodes with Degree 0: "<<numNodesD0<<endl;
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