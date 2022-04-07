//
// Created by hilto on 27/01/2021.
//

#include "Edge.h"

namespace topoSolver {

	Edge::Edge() {

	}

	Edge::Edge(int id)
	{
		//m_hed1 = hed1;
		//m_hed2 = hed2;
		m_id = id;
	}

	HalfEdge* Edge::getTwin(int hedId)
	{
		if (m_hed1->m_id == hedId)
		{
			return m_hed2;
		}
		else
		{
			return m_hed1;
		}

	}
	/* this methods just increase a already existed vector of points (recusiveley)
	*/
	void Edge::getLeafNodes(std::vector<Vertex*>* leafNodes)
	{
		if (m_isSplited)
		{
			leafNodes->push_back(m_mid);
			m_mid->markTemp = true;
			for (Edge* edgeChild : m_childrenId)
			{
				edgeChild->getLeafNodes(leafNodes);
			}
		}
	}

	void Edge::insertHed(HalfEdge* hed)
	{
		if (m_heds.empty())
		{
			m_heds.push_back(hed);
			m_hed1 = hed;
		}
		m_heds.push_back(hed);
		m_hed2 = hed;
	}

}