#include "HalfEdgeApp.h"
#include "EdgeApp.h"

namespace topoApp {
	HalfEdge::HalfEdge()
	{
	}

	HalfEdge::HalfEdge(Vertex* p1, Vertex* p2, int id)
	{
		m_p1 = p1;
		m_p2 = p2;
		m_id = id;
	}

	Vertex* HalfEdge::getP1()
	{
		return m_p1;
	}

	Vertex* HalfEdge::getP2()
	{
		return m_p2;
	}

	int HalfEdge::getId()
	{
		return m_id;
	}

	void HalfEdge::setNewPoint(int pos, Vertex* point)
	{
		if (pos == 1)
		{
			m_p1 = point;
		}
		else
		{
			m_p2 = point;
		}
	}

	void HalfEdge::initPoint(Vertex* newpt)
	{
		if (m_p1 == nullptr)
		{
			m_p1 = newpt;
			return;
		}
		m_p2 = newpt;
	}

	bool HalfEdge::isComplete()
	{
		if (!(m_p1 == nullptr) && !(m_p2 == nullptr))
		{
			return true;
		}
		return false;
	}

	void HalfEdge::setId(int id)
	{
		m_id = id;

	}

	HalfEdge* HalfEdge::getTwin()
	{
		return m_edge->getTwinHed(this);
	}



}