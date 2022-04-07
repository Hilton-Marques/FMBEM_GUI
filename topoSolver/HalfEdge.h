//
// Created by hilto on 27/01/2021.
//

#ifndef MAIN_HALFEDGE_H
#define MAIN_HALFEDGE_H


namespace topoSolver {

  class Face; // forward declaration
  class Edge; // forward declaration
  class Vertex; // forward declaration

  class HalfEdge {
  protected:
    Vertex* m_p1;
  public:

    int m_inc[2];
    int m_id = -1;
    Edge* m_edge;
    HalfEdge* m_heNext;
    Face* m_el;
    Vertex* m_p2;
    HalfEdge();
    // HalfEdge(std::vector<int> cInc, int cId, int cEdgeId, int cElId, HalfEdge* heNext);
    HalfEdge(int cInc[2], int cId, Edge* cEdge, Vertex* p1, Vertex* p2);
    Vertex* getP0() { return m_p1; };
  };
}

#endif //MAIN_HALFEDGE_H
