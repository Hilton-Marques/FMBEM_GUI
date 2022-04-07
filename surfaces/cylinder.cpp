#include "cylinder.h"
#include <cmath>
static float M_PI = 4 * atan(1);

Cylinder::Cylinder()
{
	float h_u = 2 * M_PI / m_N;
	float h_v = 1.0 / m_N;
	float teta, x, y, z;
	int count = 0;
	for (int j = 0; j < m_N; j++)
	{
		for (int i = 0; i < m_N; i++)
		{
			float tetapt[4] = { h_u * i, h_u * (i + 1), h_u * (i + 1), h_u * i };
			float zpt[4] = { (-0.5) + h_v * j, (-0.5) + h_v * j,(-0.5) + h_v * (j + 1),(-0.5) + h_v * (j + 1) };

			for (int k = 0; k < 4; k++)
			{
				teta = tetapt[k];				
				x = m_radius * cos(teta);
				y = m_radius * sin(teta);
				z = zpt[k];

				m_vertices.push_back(QVector3D(x, y, z));
				m_indices.push_back(count);
				count++;
			}

		}

	}
	std::vector<QVector3D> normals(m_vertices.size());
	for (int i = 0; i < m_vertices.size(); i++)
	{
		QVector3D p = m_vertices[i];
		normals[i] = QVector3D(p.x(), p.y(), 0.0);
	}
	m_vertices.insert(m_vertices.end(), normals.begin(), normals.end());
}
