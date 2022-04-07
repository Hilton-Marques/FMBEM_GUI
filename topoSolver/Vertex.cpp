//
// Created by hilto on 27/01/2021.
//

#include "Vertex.h"
#include "Face.h"

using namespace geomUtils;

namespace topoSolver 
{

  Vertex::Vertex() {

  }

  Vertex::Vertex(Point c_point, int c_id)
  {
    m_coord = c_point;
    m_id = c_id;
  }

  double Vertex::calculateSolidAngle()
  {
    double ret = 0.0;
    int nEl = m_elements.size();
    Face* el = m_elements[0];
    for (Face* element : m_elements)
    {
      element->markTemp = true;
    }
    for (Face* elementi : m_elements)
    {
      Point ni = elementi->getNormalUnitary();
      for (Face* elementj : elementi->m_adjacentElementsEdges)
      {
        if (elementj->markTemp) {
          Point nj = elementj->getNormalUnitary(); // normalize
          double angle = acos(ni * (-1 * nj));
          ret = ret + angle;
        }
      }
      elementi->markTemp = false;
    }
    if (nEl == 3)
      return (ret - M_PI);
    return (ret - 4 * M_PI);
  }

  void Vertex::setBC(double u)
  {
    m_bd = true;
    m_u = u;
  }

}



