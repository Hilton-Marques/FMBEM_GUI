//
// Created by hilto on 27/01/2021.
//

#include <algorithm> // find function
#include <iostream>
#include "Face.h"

namespace topoSolver {

  Face::Face() {

  }

  Face::Face(int cHedInc, double q, topoSolver::HalfEdge* heds[]) {
    m_hedInc = cHedInc;
    m_heds[0] = heds[0];
    m_heds[1] = heds[1];
    m_heds[2] = heds[2];
    if (!(q == NAN))
    {
      m_q = q;
      m_bd = true;
    }

    //getIdcPoints(heds);

  }


  std::vector<Face*> Face::getAdjacentElements(bool flag)
  {
    //std::vector<int> elIds;
    std::vector<Face*> adjacentElements;
    //adjacentElements.reserve(15);
    for (int i = 0; i < 3; i++)
    {
      HalfEdge* hed = m_heds[i];
      m_adjacentElementsEdges[i] = hed->m_edge->getTwin(hed->m_id)->m_el;
      int p2Target = hed->m_inc[1];
      // get next of twin
      hed = hed->m_edge->getTwin(hed->m_id)->m_heNext;
      while (true)
      {
        int p2 = hed->m_inc[1];
        if (p2 == p2Target)
        {
          break;
        }
        //int newId = hed->m_elId;
        topoSolver::Face* element = hed->m_el;
        if (std::find(adjacentElements.begin(), adjacentElements.end(), element) == adjacentElements.end())
        {
          adjacentElements.push_back(element);
        }
        hed = hed->m_edge->getTwin(hed->m_id)->m_heNext;
      }

    }
    //m_adjacentFaces = adjacentElements;
    return adjacentElements;
  }

  Point Face::getNormalUnitary()
  {
    return Point(m_normal / m_jacobian);
  }

  void Face::setNormal()
  {
    Point v1 = m_points[1]->m_coord - m_points[0]->m_coord;
    Point v2 = m_points[2]->m_coord - m_points[0]->m_coord;
    m_normal = Point::crossProd(v2, v1);
    m_jacobian = m_normal.getNorm();
  }

  void Face::setFMMSizes(const int& N)
  {
    m_MEG = new std::complex<double>[N];
    m_MEH = new std::complex<double>[N];
    m_N = N;

  }

  Point Face::getYc()
  {
    // Calculate centroid
    return Point((m_points[0]->m_coord + m_points[1]->m_coord + m_points[2]->m_coord) / 3.0);
  }

  void Face::getLeafElements(std::vector<Face*>* leafElements)
  {
    if (m_isLeaf)
    {
      leafElements->push_back(this);
      return;
    }
    else
    {
      for (Face* child : m_childrenElement)
      {

        child->getLeafElements(leafElements);

      }

    }


  }

  void Face::accumulate()
  {
    for (Face* child : m_childrenElement)
    {
      for (int k = 0; k < m_N; k++)
      {
        m_MEG[k] = m_MEG[k] + child->m_MEG[k];
        m_MEH[k] = m_MEH[k] + child->m_MEH[k];
      }
      child->reset();
    }
  }

  void Face::reset()
  {
    //delete[] m_MEG;
    //delete[] m_MEH;
    for (int i = 0; i < m_N; i++)
    {
      m_MEG[i] = 0.0;
      m_MEH[i] = 0.0;
    }
  }

  std::vector<Vertex*> Face::getBoundaryPoints()
  {
    std::vector<Vertex*> boundaryPoints;
    std::vector<Edge*> boundaryEdges;

    std::vector<Face*> adjacentFaces = getAdjacentElements();


    adjacentFaces.push_back(this);
    // find first Element
    HalfEdge* firstHed = nullptr;
    HalfEdge* secondHed = nullptr;
    HalfEdge* hedTwin = nullptr;
    Face* element;
    bool flag = false;
    bool mark;
    for (Face* adj : adjacentFaces)
    {
      for (HalfEdge* hed : adj->m_heds)
      {
        hedTwin = hed->m_edge->getTwin(hed->m_id);
        element = hedTwin->m_el;
        mark = !(std::find(adjacentFaces.begin(), adjacentFaces.end(), element) == adjacentFaces.end());
        if (!(mark))
        {
          firstHed = hed;
          flag = true;
          break;
        }
      }
      if (flag)
        break;
    }
    if (!flag)
    {
      return boundaryPoints;
    }
    Vertex* v0;
    Vertex* v1 = nullptr;
    v0 = firstHed->getP0();
    v0->markTemp = true;
    boundaryPoints.push_back(v0);
    boundaryEdges.push_back(firstHed->m_edge);
    while (v0 != v1)
    {
      secondHed = firstHed->m_heNext;
      element = secondHed->m_el;
      mark = !(std::find(adjacentFaces.begin(), adjacentFaces.end(), element) == adjacentFaces.end());
      while (mark)
      {
        firstHed = secondHed;
        hedTwin = firstHed->m_edge->getTwin(firstHed->m_id);
        secondHed = hedTwin->m_heNext;
        element = secondHed->m_el;
        mark = !(std::find(adjacentFaces.begin(), adjacentFaces.end(), element) == adjacentFaces.end());
      }
      v1 = firstHed->getP0();
      boundaryPoints.push_back(v1);
      v1->markTemp = true;
      boundaryEdges.push_back(firstHed->m_edge);
    }
    // get leaf nodes
    for (Edge* edge : boundaryEdges)
    {
      edge->getLeafNodes(&boundaryPoints);
    }





    return boundaryPoints;
  }

  void Face::showMEG()
  {
    std::cout << m_hedInc << "\n";
    std::cout << "MEG" << "\n";
    std::cout.precision(15);
    for (int i = 0; i < m_N; i++)
    {
      std::cout << std::imag(m_MEG[i]) << "\n";

    }
  }

  void Face::showMEH()
  {
    std::cout << m_hedInc << "\n";
    std::cout << "MEH" << "\n";
    std::cout.precision(15);
    for (int i = 0; i < m_N; i++)
    {
      std::cout << std::real(m_MEH[i]) << "\n";

    }
  }

  void Face::set_FMM_ME_sizes()
  {
    m_MEG_p.resize(3);
    m_MEH_p.resize(3);
    for (int i = 0; i < 3; i++)
    {
      m_MEH_p[i] = new std::complex<double>[m_N];
      m_MEG_p[i] = new std::complex<double>[m_N];
    }
  }

  Face::~Face()
  {
    // reset is doing the job
      //delete[] m_MEG;
      //delete[] m_MEH;
  }

}





