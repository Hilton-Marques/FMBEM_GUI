//
// Created by hilto on 27/01/2021.
//

#ifndef MAIN_FACE_H
#define MAIN_FACE_H



#include <vector>
#include <string>
#include <complex>

#include "../math/matrix.h"
#include "Vertex.h"
#include "HalfEdge.h"
#include "Edge.h"

//#include "Solid.h"

namespace topoSolver {

  class Face {

  public:
    int m_hedInc;
    int m_id;
    double m_q = NAN;
    HalfEdge* m_heds [3];
    Vertex* m_points[3];
    bool markTemp = false;
    std::complex<double>* m_MEH;

    std::complex<double>* m_MEG;

    std::vector<std::complex<double>* > m_MEG_p;

    std::vector<std::complex<double>* > m_MEH_p;


    Point m_normal;
    double m_jacobian;
    Face* m_adjacentElementsEdges[3];
    Face* m_fatherElement;
    Face* m_childrenElement[4] = { nullptr, nullptr, nullptr, nullptr };
    int m_N;
    std::vector<Vertex*> m_collectVertexes;
    bool m_isLeaf = false;
    std::vector<Face*> m_adjacentFaces;
    bool m_bd = false;

    Face();
    Face(int cHedInc, double q, topoSolver::HalfEdge* heds[]);
    std::vector<Face*> getAdjacentElements(bool flag = false);
    Point getNormalUnitary();
    void setNormal();
    void setFMMSizes(const int& N);
    Point getYc();
    void getLeafElements(std::vector<Face*>* leafElements);
    void accumulate();
    void reset();
    std::vector<Vertex*> getBoundaryPoints();
    int m_nCloseNodes;
    Matrix m_H{ 3,3 };
    Matrix m_G{ 3,3 };
    void showMEG();
    void showMEH();
    void set_FMM_ME_sizes();
    ~Face();





  };

}

#endif //MAIN_FACE_H
