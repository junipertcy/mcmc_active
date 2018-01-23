#include "Graph.h"


Graph::Graph(string sGraphFileName) {
    ifstream infile(sGraphFileName.c_str());
    string words;

    if (!infile) {
        cerr << "unable to open the input graph data file:"
             << sGraphFileName << endl;
    }
    set<unsigned int> selfloopVtxSet;
    N_ = 0;
    Q_ = 0;
    set <string> seenIds;
    set <string> seenValues;
    string word;
    string id;
    string value;
    unsigned int type;
    string sourceid;
    string targetid;
    unsigned int sourceindex;
    unsigned int targetindex;
    while (infile >> word) {
        if (word == "directed") {
            infile >> word;
            if (word == "1") {
                std::cerr << "Directed network is not supported. \n";
                throw;
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
                vtx.setIndex(N_);
                m_mapId2Index.insert(su_map_t::value_type(id, N_));
                seenIds.insert(id);
            }
            if (seenValues.count(value)) {
                type = m_mapValue2Type[value];
                vtx.setType(type);
            } else {
                vtx.setType(Q_);
                m_mapValue2Type.insert(su_map_t::value_type(value, Q_));
                seenValues.insert(value);
                Q_++;
            }
            N_++;
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
                std::cerr << "Directed network is not supported. \n";
                throw;
            }
            vtxList[sourceindex].addTarget(targetindex);
            vtxList[targetindex].addSource(sourceindex);
            vtxList[sourceindex].addSource(targetindex);
            vtxList[targetindex].addTarget(sourceindex);
        }
    }
    infile.close();

    unsigned int numEgs = 0;
    for (auto const &vtx: vtxList) {
        numEgs += vtx.getOutDegree();
    }
    E_ = numEgs / 2;
}

const Graph::Vertex &Graph::getVertex(unsigned int vtxno) const {
    if (vtxno < N_)
        return vtxList[vtxno];
    else {
        cerr << "Cannot get the vertex, the vertex no is too large." << endl
             << "-- getVertex()::Graph" << endl;
        return vtxList[N_ - 1];
    }
}

const unsigned int Graph::get_degree_at_v(unsigned int vtxno) const {
    if (vtxno > vtxList.size()) {
        cerr << "the vtx no provided is illegal.-- Graph::get_degree_at_v()" << endl;
        return 0;
    }
    return vtxList[vtxno].getInDegree();
}