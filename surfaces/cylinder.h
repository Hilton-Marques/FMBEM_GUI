#ifndef CYLINDER_H
#define CYLINDER_H

#include <vector>
#include <QVector3D>
#include <QOpenGLFunctions>

class Cylinder
{
public:
	Cylinder();
	const int m_N = 20;
	qreal m_radius = 0.012;

	std::vector<QVector3D> m_vertices;
	std::vector<GLushort> m_indices;

private:

};
#endif