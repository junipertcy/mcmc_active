#ifndef GRAPH_H
#define GRAPH_H

#include "Global.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <set>

class Graph {
private:
    class Vertex;

public:
    explicit Graph(string sGraphFile);

    uint_vec_t vtxSelfloopFlag;

    unsigned getNumVtx() const { return m_numVtx; }

    unsigned getNumType() const { return m_numType; }

    void setDirected(bool b_direc) { m_directed = b_direc; }

    void setSelfloop(bool b_sl) { m_selfloop = b_sl; }

    bool isDirected() const { return m_directed; }

    bool hasSelfloop() const { return m_selfloop; }

    const Vertex &getVertex(unsigned vtxno) const;

    bool vtxHasSelfloop(unsigned vtxno) const { return vtxSelfloopFlag[vtxno]; }

    const unsigned getVtxDegree(unsigned vtxno) const;

private:
    void calcNumEdges();

    void calcGroupConnNumMatrix();

    void calcGroupConnMatrix();

    unsigned m_numVtx;
    unsigned m_numType;
    unsigned m_numEgs;
    bool m_directed;
    bool m_selfloop;
    unsigned m_numSelfloop;
    vector <Vertex> vtxList;
    map<unsigned, string> m_mapIndex2Id;
    map<string, unsigned> m_mapId2Index;
    map<string, unsigned> m_mapValue2Type;
    map<unsigned, string> m_mapType2Value;

    uint_mat_t m_groupConnNumMatrix;
    float_mat_t m_groupConnMatrix;
    uint_vec_t m_groupCardi;

    class Vertex {
        friend class Graph;

    public:
        Vertex() {
            m_id="";
            m_value="";
            m_index=0;
            m_type=0;
        }
        string m_id;
        string m_value;
        unsigned int m_index;
        unsigned int m_type;

        set<unsigned int> m_targets;
        set<unsigned int> m_sources;

        void setId(string vtx_id) { m_id = vtx_id; }

        void setIndex(unsigned int vtx_index) { m_index = vtx_index; }

        void setValue(string vtx_value) { m_value = vtx_value; }

        void setType(unsigned int vtx_type) { m_type = vtx_type; }

        unsigned getType() const { return m_type; }

        void addTarget(unsigned int target) { m_targets.insert(target); }

        void addSource(unsigned int source) { m_sources.insert(source); }

        const set<unsigned int> &getSources() const { return m_sources; }

        const set<unsigned int> &getTargets() const { return m_targets; }

        const unsigned int getOutDegree() const { return m_targets.size(); }

        const unsigned int getInDegree() const { return m_sources.size(); }
    };
};

#endif
