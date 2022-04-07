#include "cone.h"
#include <cmath>
static const double M_PI = 4 * std::atan(1);
Cone::Cone()
{
	float h_u = 2 * M_PI / m_N;
	float h_v = (m_radius) / m_N;
	float teta, x, y, z, alpha,t;
	alpha = -1.0 / m_radius;
	int count = 0;
	for (int j = 0; j < m_N; j++)
	{
		for (int i = 0; i < m_N; i++)
		{
			float tetapt[4] = { h_u * i, h_u * (i + 1), h_u * (i + 1), h_u * i };
			float tpt[4] = { h_v * j, h_v * j, h_v * (j + 1),h_v * (j + 1) };

			for (int k = 0; k < 4; k++)
			{
				teta = tetapt[k];
				t = tpt[k];
				x = t * -sin(teta);
				y = t * cos(teta);
				z = (0.5 + alpha*t);

				m_vertices.push_back(QVector3D(x, y, z));
				m_indices.push_back(count);
				count++;
			}

		}
	}
}
