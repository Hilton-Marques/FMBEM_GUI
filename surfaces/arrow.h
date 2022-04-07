#ifndef ARROW_H
#define ARROW_H
#include <vector>
#include <QVector3D>
#include <QOpenGLFunctions>
class Arrow
{
public:
	Arrow();
	std::vector<QVector3D> m_vertices;
	std::vector<GLushort> m_indices;
};



#endif