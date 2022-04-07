//
// Created by hilto on 27/01/2021.
//

#include "HalfEdge.h"
namespace topoSolver {

  HalfEdge::HalfEdge() {

  }
  HalfEdge::HalfEdge(int cInc[2], int cId, Edge* cEdge, Vertex* p1, Vertex* p2)
  {
    m_inc[0] = cInc[0], m_inc[1] = cInc[1];
    m_id = cId;
    m_edge = cEdge;
    m_p1 = p1, m_p2 = p2;
  }
}




