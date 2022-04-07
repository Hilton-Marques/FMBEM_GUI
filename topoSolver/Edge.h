//
// Created by hilto on 27/01/2021.
//

#ifndef MAIN_EDGE_H
#define MAIN_EDGE_H

#include "HalfEdge.h"
#include "Vertex.h"

namespace topoSolver
{

  class Edge {
  public:
    bool m_isSplited = false;
    HalfEdge* m_hed1;
    HalfEdge* m_hed2;
    Edge* m_childrenId[2];
    int m_id;
    Vertex* m_mid;
    std::vector<HalfEdge*> m_heds;
    Edge();
    Edge(int id);
    HalfEdge* getTwin(int hedId);
    void getLeafNodes(std::vector<Vertex*>* leafNodes);
    void insertHed(HalfEdge* hed);

  };

}

#endif //MAIN_EDGE_H
