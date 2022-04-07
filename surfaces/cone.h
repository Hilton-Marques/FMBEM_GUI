#ifndef CONE_H
#define CONE_H
#include <vector>
#include <QVector3D>
#include <QOpenGLFunctions>
class Cone
{
public:
	Cone();
	const int m_N = 20;
	const qreal m_radius = 1;
	const qreal m_tronco = 0.2;
	std::vector<QVector3D> m_vertices;
	std::vector<GLushort> m_indices;

private:

};




#endif // !CONE_H

