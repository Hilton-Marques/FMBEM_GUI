#include "convexHull.h"

ConvexHull::ConvexHull(Model* model)
{
  m_model = model;
  m_verts = model->m_vertexes;
}

void ConvexHull::initModel()
{
  std::vector<HalfEdge*> tempHeds;
  Face face1 = findTriangleOnHull();
  for (int i = 1; i < 3; i++)
  {
    HalfEdge* hedTwin = face1.get3Heds()[i]->getTwin();
    tempHeds.push_back(hedTwin);
  }
  while (!tempHeds.empty())
  {
    HalfEdge* hedi = tempHeds.back();
    tempHeds.pop_back();
    Edge* edgei = hedi->getEdge();
    
    if (!edgei->isProcessed())
    {
      Vertex* q = pivotOnEdge(hedi);
      int incs[2][2] = { {hedi->getP2()->getId(),q->getId()}, {q->getId(), hedi->getP1()->getId()} };
      std::vector<HalfEdge*> hedsTri;
      for (int i = 0; i < 2; i++)
      {
        bool flag = false;
        #pragma omp parallel for
        for (int j = 0; j < m_model->m_edges.size(); j++)
        {
          Edge* edge = m_model->m_edges[j];
          if (incs[i][0] == edge->getIdP1() && incs[i][1] == edge->getIdP2() ||
            incs[i][1] == edge->getIdP1() && incs[i][0] == edge->getIdP2())
          {
            flag = true;
            edge->setIsProcessed(true);
            if (incs[i][0] == edge->getIdP1())
            {
              hedsTri.push_back(edge->getCurrentHed());
            }
            else
            {
              hedsTri.push_back(edge->getTwinHed());
            }              
            break;
          }
        }
        if (!flag)
        {
          // create new Edges
            //crate new Faces
           // edge1
          HalfEdge* hed = new HalfEdge(m_model->m_vertexes[incs[i][0]], m_model->m_vertexes[incs[i][1]], m_model->m_heds.size());
          m_model->m_heds.push_back(hed);
          HalfEdge* hedTwin = new HalfEdge(m_model->m_vertexes[incs[i][1]], m_model->m_vertexes[incs[i][0]], m_model->m_heds.size());
          m_model->m_heds.push_back(hedTwin);
          Edge* edge = new Edge(hed, hedTwin, m_model->m_edges.size());
          m_model->m_edges.push_back(edge);
          hed->setEdge(edge);
          hedTwin->setEdge(edge);
          tempHeds.push_back(hedTwin);
          hedsTri.push_back(hed);
        }
      }
      Face* t = new Face(hedi, hedsTri[0], hedsTri[1], m_model->m_faces.size());
      m_model->m_faces.push_back(t);
      edgei->setIsProcessed(true);
    }
  }
}

Face ConvexHull::findTriangleOnHull()
{
  HalfEdge* hed0 = findEdgeOnHull();
  Vertex* r = pivotOnEdge(hed0);
  //crate new Faces
  // edge1
  HalfEdge* hed1 = new HalfEdge(hed0->getP2(), r, m_model->m_heds.size());
  m_model->m_heds.push_back(hed1);
  HalfEdge* hedTwin = new HalfEdge(r, hed0->getP2(), m_model->m_heds.size());
  m_model->m_heds.push_back(hedTwin);
  Edge* edge = new Edge(hed1, hedTwin, m_model->m_edges.size());
  m_model->m_edges.push_back(edge);
  hed1->setEdge(edge);
  hedTwin->setEdge(edge);
  // edge2
  HalfEdge* hed2 = new HalfEdge(r, hed0->getP1(), m_model->m_heds.size());
  m_model->m_heds.push_back(hed2);
  hedTwin = new HalfEdge(hed0->getP1(), r, m_model->m_heds.size());
  m_model->m_heds.push_back(hedTwin);
  edge = new Edge(hed2, hedTwin, m_model->m_edges.size());
  m_model->m_edges.push_back(edge);
  hed2->setEdge(edge);
  hedTwin->setEdge(edge);
  // new Face
  Face t = Face(hed0, hed1, hed2, 42);
  m_model->m_faces.push_back( new Face(hed0, hed1, hed2, m_model->m_faces.size()) );
  
  return t;
}

Vertex* ConvexHull::pivotOnEdge(HalfEdge* hed)
{
  
  QVector3D* q0 = hed->getP1()->getP();
  QVector3D* q1 = hed->getP2()->getP();
  Vertex* r = m_verts[0];
  QVector3D* p = r->getP();
  qreal area2 = squaredArea(q0, q1, p);
  for (int i = 1; i < m_verts.size(); i++)
  {
    QVector3D* pi = m_verts[i]->getP();
    qreal volume = signedVolume(q0, q1, p, pi);
    if (volume < 0)
    {
      p = pi;
      r = m_verts[i];
    }
    else if( volume == 0) {
      qreal area2_ = squaredArea(q0, q1, pi);
      if (area2_ > area2)
      {
        p = pi;
        r = m_verts[i];
        area2 = area2_;
      }

    }
  }
  return r;
}

HalfEdge* ConvexHull::findEdgeOnHull()
{
  //get most x point
  Vertex* ptRight = m_verts[0];
  qreal xRight = ptRight->getP()->x();
  for (int i = 1; i < m_verts.size(); i++)
  {
    auto pt = m_verts[i];
    qreal xRight_ = pt->getP()->x();
    if (xRight_ > xRight )
    {
      xRight = xRight_;
      ptRight = pt;
    }
  }
  QVector3D supportCoord = QVector3D(ptRight->getP()->x() , ptRight->getP()->y() + 1.0, ptRight->getP()->z());
  Vertex supportV(supportCoord, m_verts.size());
  HalfEdge supportHed (ptRight, &supportV , 0);
  Vertex* r = pivotOnEdge(&supportHed);
  //create edge
  HalfEdge* hed1 = new HalfEdge(ptRight, r, m_model->m_heds.size());
  m_model->m_heds.push_back(hed1);
  HalfEdge* hed2 = new HalfEdge(r, ptRight, m_model->m_heds.size());
  m_model->m_heds.push_back(hed2);
  Edge* edge = new Edge(hed1, hed2, m_model->m_edges.size());
  m_model->m_edges.push_back(edge);
  hed1->setEdge(edge);
  hed2->setEdge(edge);
  return hed1;
}

qreal ConvexHull::squaredArea(QVector3D* p0, QVector3D* p1, QVector3D* p2)
{
  QVector3D u = *p1 - *p0;
  QVector3D v = *p2 - *p0;
  QVector3D area = QVector3D::crossProduct(u, v);
  return QVector3D::dotProduct(area,area);
}

qreal ConvexHull::signedVolume(QVector3D* p0, QVector3D* p1, QVector3D* p2, QVector3D* p3)
{
  QVector3D u = *p1 - *p0;
  QVector3D v = *p2 - *p0;
  QVector3D z = *p3 - *p0;
  float data[] = {   u.x(),u.y(),u.z() ,0,
                     v.x(),v.y(),v.z(),  0,
                     z.x(),z.y(),z.z(), 0,
                     0.0,0.0,0.0,1.0  };
  QMatrix4x4 matrix (data);
  float det = matrix.determinant();
  return std::abs(det) < 0.00000001 ? 0 : det;
}
