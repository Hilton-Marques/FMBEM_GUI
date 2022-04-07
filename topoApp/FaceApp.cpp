#include "FaceApp.h"

namespace topoApp {
	Face::Face()
	{
	}

	Face::Face(HalfEdge* e1, HalfEdge* e2, HalfEdge* e3, int id)
	{
		m_he1 = e1;
		m_he2 = e2;
		m_he3 = e3;
		m_id = id;

		m_p1 = m_he1->getP1()->getP();
		m_p2 = m_he1->getP2()->getP();
		m_p3 = m_he2->getP2()->getP();

		m_color = QVector3D(0.0, 0.0, 1.0);







	}

	QVector3D Face::getCentroide()
	{
		return QVector3D((*m_p1 + *m_p2 + *m_p3) / 3);
	}

	std::vector<QVector3D> Face::getPoints()
	{
		std::vector<QVector3D> out = { *m_p1, *m_p2, *m_p3 };
		return out;

	}

	bool Face::isScreenPointClose(QVector2D point, QMatrix4x4* transform)
	{
		QVector3D centroide = getCentroide();
		//get depth
		QVector4D centroideClipping = (*transform) * QVector4D(centroide, 1.0);
		qreal a = centroideClipping.x() / centroideClipping.w();
		qreal b = centroideClipping.y() / centroideClipping.w();

		m_depth = centroideClipping.z() / centroideClipping.w();
		// start
		std::vector<QVector3D*> points = { m_p1, m_p2, m_p3 };
		std::vector<QVector4D> pointsClipping(3);
		std::vector<QVector2D> pointsClippingNDC(3);
		for (int i = 0; i < 3; i++)
		{
			pointsClipping[i] = (*transform) * QVector4D(*points[i], 1.0);
			pointsClippingNDC[i] = QVector2D(pointsClipping[i].x() / pointsClipping[i].w(), pointsClipping[i].y() / pointsClipping[i].w());
		}
		// check if order was reversed with camera
		if (orientedArea(pointsClippingNDC[0], pointsClippingNDC[1], pointsClippingNDC[2]) < 0)
		{
			pointsClippingNDC = { pointsClippingNDC[0] , pointsClippingNDC[2], pointsClippingNDC[1] };
		}
		QVector2D p1, p2, v, u;
		for (int i = 0; i < 3; i++)
		{
			p1 = pointsClippingNDC[i];
			p2 = pointsClippingNDC[(i + 1) % 3];
			if (orientedArea(p1, p2, point) < 0)
			{
				return false;
			}
		}
		//reverse order

		return true;
	}

	QVector3D Face::getNormalColor()
	{
		return QVector3D(0.0, 0.0, 1.0); //blue
	}

	QVector3D Face::getSelectedColor()
	{
		return QVector3D(0.0, 0.4, 0.4);
	}

	QVector3D Face::getCurrentColor()
	{
		return m_color;
	}

	void Face::setColor(QVector3D color)
	{
		m_color = color;
	}

	qreal Face::getDepth()
	{
		return m_depth;
	}

	qreal Face::orientedArea(QVector2D a, QVector2D b, QVector2D c)
	{
		QVector2D v = (b - a).normalized();
		QVector2D u = (c - a).normalized();
		QVector2D vPerp(-v.y(), v.x());

		return QVector2D::dotProduct(u, vPerp);
	}




	int Face::getId()
	{
		return m_id;
	}

	void Face::setId(int id)
	{
		m_id = id;
	}

	Edge* Face::getE1()
	{
		return m_he1->getEdge();
	}

	Edge* Face::getE2()
	{
		return m_he2->getEdge();
	}

	Edge* Face::getE3()
	{
		return m_he3->getEdge();
	}

	HalfEdge* Face::getValidHed(Edge* edge)
	{
		if (m_heds.empty())
		{
			return edge->getCurrentHed();
		}
		HalfEdge* lastHed = m_heds.back();
		if (edge->getCurrentHed()->getP1() == lastHed->getP2())
		{
			return edge->getCurrentHed();
		}
		else if (edge->getTwinHed()->getP1() == lastHed->getP2())
		{
			return edge->getTwinHed();
		}
		return nullptr;
	}

	void Face::insertHed(HalfEdge* hed)
	{
		m_heds.push_back(hed);
	}

	bool Face::isComplete()
	{
		int n = m_heds.size();
		if (n < 3)
			return false;
		HalfEdge* firstHed = m_heds[0];
		HalfEdge* lastHed = m_heds.back();
		if (firstHed->getP1() == lastHed->getP2())
		{
			return true;
		}
		return false;
	}

	void Face::reset()
	{
		for (HalfEdge* hed : m_heds)
		{
			hed->getEdge()->setColorEdge(Edge::normalColor);
		}
	}

	bool Face::contain(Edge* edge)
	{
		if (getE1() == edge || getE2() == edge || getE3() == edge)
		{
			return true;
		}

		return false;
	}

	QVector3D Face::getNormal()
	{
		QVector3D normal = -QVector3D::crossProduct(*m_p2 - *m_p1, *m_p3 - *m_p1);
		return normal.normalized();
	}

	void Face::reverseHeds()
	{
		int n = m_heds.size();
		std::vector<HalfEdge*> newHeds(n);
		int count = 0;
		for (int i = n - 1; i > -1; i--)
		{
			HalfEdge* hed = m_heds[i];
			newHeds[count] = hed->getEdge()->getTwinHed(hed);
			count++;
		}
		m_heds = newHeds;
	}

	qreal Face::getArea()
	{
		QVector3D normal = QVector3D::crossProduct(*m_p2 - *m_p1, *m_p3 - *m_p1);
		return 0.5 * normal.length();
	}

	std::vector<HalfEdge*> Face::get3Heds()
	{
		std::vector<HalfEdge*> out = { m_he1, m_he2, m_he3 };
		return out;
	}

	std::vector<QVector3D> Face::getPointsAndNormals()
	{
		std::vector<QVector3D> out(6);
		QVector3D n = getNormal();
		out[0] = *m_p1;
		out[1] = *m_p2;
		out[2] = *m_p3;
		out[3] = n;
		out[4] = n;
		out[5] = n;
		return out;
	}
}


