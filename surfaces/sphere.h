#ifndef SPHERE_H
#define SPHERE_H
#include <vector>
#include <QOpenGLFunctions>
#include <QVector3D>
class Sphere
{
public:
	Sphere();
	int m_N = 20 ;
	qreal m_radius = 0.055;
	std::vector<QVector3D> m_vertices;
	std::vector<GLushort> m_indices;

private:

};



#endif // !SPHERE_H
