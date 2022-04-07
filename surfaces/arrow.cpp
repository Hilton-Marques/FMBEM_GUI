#include "arrow.h"
#include "cone.h"
#include "cylinder.h"
#include <QVector2D>
Arrow::Arrow()
{
  Cone cone;
  Cylinder cylinder;
  int count = 0;
  for (int i = 0; i < cylinder.m_indices.size(); i++)
  {
    QVector3D pt = cylinder.m_vertices[i];
    qreal x = pt.x() ;
    qreal y = pt.y() ;
    qreal z = 0.25 * (pt.z() + 0.5) ;
    m_vertices.push_back(QVector3D(x, y, z));
    m_indices.push_back(count);
    count++;
  }
  for (int i = 0; i < cone.m_indices.size(); i++)
  {
    QVector3D pt = cone.m_vertices[i];
    qreal x = pt.x() * 0.025;
    qreal y = pt.y() * 0.025;
    qreal z = 0.1* pt.z() + 0.25;
    m_vertices.push_back(QVector3D(x,y,z));
    m_indices.push_back(count);
    count++;
  }
  std::vector<QVector3D> normals(m_vertices.size());
  for (int i = 0; i < cylinder.m_indices.size(); i++)
  {
    QVector3D pt = m_vertices[i];
    normals[i] = QVector3D(pt.x(), pt.y(), 0.0);
  }
  for (int i = cylinder.m_indices.size(); i < m_vertices.size(); i++)
  {
    QVector3D pt = m_vertices[i];
    QVector2D proj = QVector2D ( pt.x(), pt.y() );
    normals[i] = QVector3D(pt.x(), pt.y(), proj.length());
  }
  m_vertices.insert(m_vertices.end(), normals.begin(), normals.end());

}
