#include "EdgeApp.h"
#include <QVector4D>
#include <QVector2D>

namespace topoApp {
  const QVector3D Edge::normalColor = QVector3D(0.0, 1.0, 0.0);
  const QVector3D Edge::selectedColor = QVector3D(0.5, 0.5, 0.5);

  Edge::Edge()
  {
  }

  Edge::Edge(HalfEdge* hed, HalfEdge* twin, int id)
  {
    m_hed = hed;
    m_twinHed = twin;
    m_centroide = 0.5 * (*p1() + *p2());
    m_colorEdge = QVector3D(0.0, 1.0, 0.0);
    m_colorEdge = QVector3D(0.0, 1.0, 0.4);
    m_id = id;
    m_p1 = m_hed->getP1();
    m_p2 = m_hed->getP2();
    createLocalCoordinateSystem();
    buildBB();
  }

  QVector3D* Edge::p1()
  {
    return m_hed->getP1()->getP();
  }

  QVector3D* Edge::p2()
  {
    return m_hed->getP2()->getP();
  }

  void Edge::buildBB()
  {
    m_bb.resize(8);
    QVector3D edge;
    edge = (*p2() - *p1());
    m_bb[0] = *p1() - 1.1 * m_radius * m_x; //left
    m_bb[1] = *p1() + 1.1 * m_radius * m_x; //left
    m_bb[2] = *p1() - 1.1 * m_radius * m_y; //left
    m_bb[3] = *p1() + 1.1 * m_radius * m_y; //left

    m_bb[4] = m_bb[0] + edge; //left
    m_bb[5] = m_bb[1] + edge; //left
    m_bb[6] = m_bb[2] + edge; //left
    m_bb[7] = m_bb[3] + edge; //left


  }

  void Edge::getBBOnClippingCoordinates(qreal* xmin, qreal* xmax, qreal* ymin, qreal* ymax, QMatrix4x4* transform)
  {
    //get depth of centroide
    QVector4D centroideClipping = (*transform) * QVector4D(m_centroide, 1.0);
    m_depth = centroideClipping.z() / centroideClipping.w();

    QVector3D p;
    QVector4D pClipping;
    qreal xS, yS;
    qreal xminP0, xmaxP0, yminP0, ymaxP0;
    qreal xminP1, xmaxP1, yminP1, ymaxP1;
    p = m_bb[0];
    pClipping = (*transform) * QVector4D(p, 1.0);
    xS = pClipping.x() / pClipping.w();
    yS = pClipping.y() / pClipping.w();
    xminP0 = xS, xmaxP0 = xS, yminP0 = yS, ymaxP0 = yS;
    // On this bb we have the minimum point on left and the max on right
    for (size_t i = 1; i < 4; i++)
    {
      p = m_bb[i];
      pClipping = (*transform) * QVector4D(p, 1.0);
      xS = pClipping.x() / pClipping.w();
      yS = pClipping.y() / pClipping.w();
      if (xS < *xmin)
      {
        *xmin = xS;
      }
      else if (xS > *xmax)
      {
        *xmax = xS;
      }
      if (yS < *ymin)
      {
        *ymin = yS;
      }
      else if (yS > *ymax)
      {
        *ymax = yS;
      }
    }




  }

  qreal Edge::getDepth()
  {
    return m_depth;
  }

  bool Edge::isScreenPointClose(QVector2D point, QMatrix4x4* transform)
  {
    //get depth
    QVector4D centroideClipping = (*transform) * QVector4D(m_centroide, 1.0);
    m_depth = centroideClipping.z() / centroideClipping.w();
    // start
    QVector4D p1Clipping, p2Clipping;
    p1Clipping = (*transform) * QVector4D(*p1(), 1.0);
    p2Clipping = (*transform) * QVector4D(*p2(), 1.0);
    QVector2D p1ClippingNDC(p1Clipping.x() / p1Clipping.w(), p1Clipping.y() / p1Clipping.w());
    QVector2D p2ClippingNDC(p2Clipping.x() / p2Clipping.w(), p2Clipping.y() / p2Clipping.w());
    QVector2D v = (p2ClippingNDC - p1ClippingNDC);
    QVector2D vUnit = v.normalized();
    QVector2D u = point - p1ClippingNDC;
    qreal projParallelLength = (QVector2D::dotProduct(u, vUnit));
    if (projParallelLength > v.length() || projParallelLength < 0)
    {
      return false;
    }
    QVector2D projParallel = projParallelLength * vUnit;
    QVector2D projOrtho = u - projParallel;
    qreal d = projOrtho.lengthSquared();

    // find screen tol
    QVector3D p;
    QVector4D pClipping;
    QVector2D pClippingNDC;
    qreal screenTol, screenTolTemp;
    p = m_bb[0];
    pClipping = (*transform) * QVector4D(p, 1.0);
    pClippingNDC = QVector2D(pClipping.x() / pClipping.w(), pClipping.y() / pClipping.w());
    screenTol = (pClippingNDC - p1ClippingNDC).lengthSquared();
    for (int i = 1; i < 4; i++)
    {
      p = m_bb[i];
      pClipping = (*transform) * QVector4D(p, 1.0);
      pClippingNDC = QVector2D(pClipping.x() / pClipping.w(), pClipping.y() / pClipping.w());
      screenTolTemp = (pClippingNDC - p1ClippingNDC).lengthSquared();
      if (screenTolTemp > screenTol)
      {
        screenTol = screenTolTemp;
      }
    }

    return (d <= screenTol);

  }

  QVector3D Edge::getColorEdge()
  {
    return m_colorEdge;
  }

  QVector3D Edge::getColorCone()
  {
    return m_colorEdge;
  }

  QVector3D Edge::getNormalColorEdge()
  {
    return QVector3D(0.0, 1.0, 0.0); // green
  }

  QVector3D Edge::getNormalColorCone()
  {
    return QVector3D(0.0, 1.0, 0.4); // yellow
  }

  QVector3D Edge::getSelectedColorEdge()
  {
    return QVector3D(0.5, 0.5, 0.0);
  }

  QVector3D Edge::getSelectedColorCone()
  {
    return QVector3D(0.3, 0.3, 0.0);
  }

  void Edge::setColorEdge(QVector3D color)
  {
    m_colorEdge = color;
  }

  void Edge::setColorCone(QVector3D color)
  {
    m_colorCone = color;
  }

  void Edge::createLocalCoordinateSystem()
  {
    QVector3D edge = *p2() - *p1();
    m_z = edge.normalized();
    m_x = QVector3D(-m_z.y(), m_z.x(), m_z.z());
    if (QVector3D::dotProduct(m_x, m_z) == 1)
    {
      m_x = QVector3D(m_z.x(), m_z.z(), -m_z.y());
    }
    m_y = QVector3D::crossProduct(m_z, m_x);
  }

  int Edge::getId()
  {
    return m_id;
  }

  int Edge::getIdP1()
  {
    return m_hed->getP1()->getId();
  }

  int Edge::getIdP2()
  {
    return m_hed->getP2()->getId();
  }

  void Edge::setNewPoint(int pos, Vertex* point)
  {
    m_hed->setNewPoint(pos, point);
    // update
    createLocalCoordinateSystem();
    buildBB();
  }

  HalfEdge* Edge::getCurrentHed()
  {
    return m_hed;
  }

  HalfEdge* Edge::getTwinHed()
  {
    return m_twinHed;
  }

  HalfEdge* Edge::getTwinHed(HalfEdge* hed)
  {
    if (hed == m_hed)
    {
      return m_twinHed;
    }
    else
    {
      return m_hed;
    }
  }

  bool Edge::contain(Vertex* v)
  {
    if ((m_p1 == v) || (m_p2 == v))
    {
      return true;
    }
    return false;
  }
}