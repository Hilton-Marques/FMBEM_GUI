//
// Created by hilto on 27/01/2021.
//

#ifndef MAIN_VERTEX_H
#define MAIN_VERTEX_H

#include "../geom/Point.h"
#include <vector>

namespace topoSolver {
  class Face;
  class Vertex {

  public:
    Vertex();
    Vertex(Point c_point, int c_id);
    double calculateSolidAngle();
    void setBC(double u); // set boundary condition
    Point m_coord;
    int m_id;
    double m_u = NAN;
    bool markTemp = false;
    bool m_bd = false;
    std::vector <Face*> m_elements;
  };

}

#endif //MAIN_VERTEX_H
