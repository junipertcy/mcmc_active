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

    unsigned int get_N() const { return N_; }

    unsigned int get_Q() const { return Q_; }

    const Vertex &getVertex(unsigned vtxno) const;

    const unsigned int get_degree_at_v(unsigned vtxno) const;

    ~Graph();

private:
    unsigned int N_;
    unsigned int Q_;
    unsigned int E_;

    vector <Vertex> vtxList;
    su_map_t m_mapId2Index;
    su_map_t m_mapValue2Type;

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
