#include "sphere.h"
#include <cmath>
static float M_PI = 4 * atan(1);
Sphere::Sphere()
{
	std::vector<QVector3D> normals;
	float h_u = 2 * M_PI / m_N;
	float h_v = M_PI / m_N;
	float fi, teta ,x,y,z;
	int count = 0;
	for (int j = 0; j < m_N; j++)
	{
		for (int i = 0; i < m_N; i++)
		{
			float tetapt[4] = { h_u * i, h_u * (i + 1), h_u * (i + 1), h_u * i };
			float fipt[4] = { (-M_PI / 2) + h_v * j, (-M_PI / 2)  + h_v * j,(-M_PI / 2) + h_v * (j + 1),(-M_PI / 2) + h_v * (j + 1) };
			//int ids[4] = { count, count+1,count+2,count+3 };
			//count = count + 4;
			QVector3D points[4];
			for (int k = 0; k < 4; k++)
			{
				teta = tetapt[k];
				fi =  fipt[k];
				x = m_radius * cos(teta) * cos(fi);
				y = m_radius * sin(teta) * cos(fi);
				z = m_radius * sin(fi);
				QVector3D pt = QVector3D(x, y, z);
				m_vertices.push_back(pt );
				m_indices.push_back(count);
				count++;
				points[k] = pt;
			}
			QVector3D normal = QVector3D::crossProduct(points[1] - points[0], points[2] - points[0]);
			normals.push_back(normal);
		}
	}
	m_vertices.insert(m_vertices.end(), m_vertices.begin(), m_vertices.end());
}
