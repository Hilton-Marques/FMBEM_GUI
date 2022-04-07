
#include "VertexApp.h"

//fixed color

namespace topoApp {
  const QVector3D Vertex::normalColor = QVector3D(1.0, 0.0, 0.0);
  const QVector3D Vertex::selectedColor = QVector3D(1.0, 0.0, 1.0);

  Vertex::Vertex()
  {
  }

  Vertex::Vertex(QVector3D& p, int id)
  {
    m_p = p;
    m_color = QVector3D(1.0, 0.0, 0.0); // yellow
    m_id = id;
    buildBB();
  }

  QVector3D* Vertex::getP()
  {
    return &m_p;
  }

  QVector3D Vertex::getColor()
  {
    return m_color;
  }

  QVector3D Vertex::getSelectedColor()
  {
    return QVector3D(1.0, 0.0, 1.0);
  }

  QVector3D Vertex::getNormalColor()
  {
    return QVector3D(1.0, 0.0, 0.0);
  }

  qreal Vertex::getDepth()
  {
    return m_depth;
  }

  std::vector<QVector3D>* Vertex::getBB()
  {
    return &m_bb;
  }

  void Vertex::buildBB()
  {

    QVector3D left, right, bottom, top, front, back;
    m_bb.resize(6);
    m_bb[0] = left = m_radius * QVector3D(-1.0, 0.0, 0.0) + m_p;
    m_bb[1] = right = m_radius * QVector3D(1.0, 0.0, 0.0) + m_p;
    m_bb[2] = bottom = m_radius * QVector3D(0.0, -1.0, 0.0) + m_p;
    m_bb[3] = top = m_radius * QVector3D(0.0, 1.0, 0.0) + m_p;
    m_bb[4] = front = m_radius * QVector3D(0.0, 0.0, 1.0) + m_p;
    m_bb[5] = back = m_radius * QVector3D(0.0, 0.0, -1.0) + m_p;

  }
  void Vertex::setColor(QVector3D color)
  {
    m_color = color;
  }
  void Vertex::getBBOnClippingCoordinates(qreal* xmin, qreal* xmax, qreal* ymin, qreal* ymax, QMatrix4x4* transform)
  {

    QVector3D p;
    QVector4D pClipping;
    qreal xS, yS;
    p = m_bb[0];
    pClipping = (*transform) * QVector4D(p, 1.0);
    xS = pClipping.x() / pClipping.w();
    yS = pClipping.y() / pClipping.w();
    m_depth = pClipping.z() / pClipping.w();
    *xmin = xS, * xmax = xS, * ymin = yS, * ymax = yS;


    for (size_t i = 1; i < m_bb.size(); i++)
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

  bool Vertex::isScreenPointClose(const QVector2D point, QMatrix4x4* transform)
  {
    QVector4D p1Clipping;
    //qreal screenTol = 0.0005;
    p1Clipping = (*transform) * QVector4D(m_p, 1.0);
    m_depth = p1Clipping.z() / p1Clipping.w();
    QVector2D p1ClippingNDC(p1Clipping.x() / p1Clipping.w(), p1Clipping.y() / p1Clipping.w());
    QVector2D v = point - p1ClippingNDC;

    // find screen tol
    QVector3D p;
    QVector4D pClipping;
    QVector2D pClippingNDC;
    qreal screenTol, screenTolTemp;
    p = m_bb[0];
    pClipping = (*transform) * QVector4D(p, 1.0);
    pClippingNDC = QVector2D(pClipping.x() / pClipping.w(), pClipping.y() / pClipping.w());
    screenTol = (pClippingNDC - p1ClippingNDC).lengthSquared();
    for (int i = 1; i < m_bb.size(); i++)
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

    return (v.lengthSquared() <= screenTol);
  }

  int Vertex::getId()
  {
    return m_id;
  }

  qreal Vertex::X()
  {
    return m_p.x();
  }

  qreal Vertex::Y()
  {
    return m_p.y();
  }

  qreal Vertex::Z()
  {
    return m_p.z();
  }

  void Vertex::setPosValue(int pos, qreal value)
  {
    switch (pos)
    {
    case 0:
      m_p.setX(value);
      break;
    case 1:
      m_p.setY(value);
      break;
    case 2:
      m_p.setZ(value);
      break;
    }
    // update bounding box
    buildBB();
  }
}
