#include "collector.h"

Collector::Collector()
{
}

Collector::~Collector()
{
	delete m_hedBeingCollected;
	delete m_faceBeingCollected;
}

void Collector::startCollection()
{
	m_startCollecting = true;
	switch (m_objectType)
	{
	case ObjectType::UNDEF:
		break;
	case ObjectType::HALFEDGE:
		delete m_hedBeingCollected;
		m_hedBeingCollected = new HalfEdge();
		break;
	case ObjectType::FACE:
		delete m_faceBeingCollected;
		m_faceBeingCollected = new Face();
		break;
	}
}

void Collector::setObjectType(ObjectType type)
{
	switch (type)
	{
	case ObjectType::UNDEF:
		break;
	case ObjectType::HALFEDGE:
		m_objectType = type;
		break;
	case ObjectType::FACE:
		m_objectType = type;
		break;
	}
}

ObjectType Collector::getObjectType()
{
	return m_objectType;
}

bool Collector::isPointSelected(QVector2D &pt, QMatrix4x4* transform,const Model* model)
{
		qreal screenTol = 0.0006;
		qreal minDepth = 1.0;
		int idSelected = -1;

		for (int i = 0; i < model->m_vertexes.size(); i++)
		{
			Vertex* v = model->m_vertexes[i];
			if (v->isScreenPointClose(pt, transform))
			{
				if (v->getDepth() < minDepth)
				{
					minDepth = v->getDepth();
					idSelected = i;
				}
			}
		}
		if (!(idSelected == -1))
		{
			Vertex* newVertex = model->m_vertexes[idSelected];
			newVertex->setColor(Vertex::selectedColor);
			m_hedBeingCollected->initPoint(newVertex);
			return true;
		}
		return false;
}

bool Collector::isEdgeSelected(QVector2D& pt, QMatrix4x4* transform, const Model* model)
{
	int idSelected = -1;
	qreal minDepth = 1.0;

	// Edges
	for (int i = 0; i < model->m_edges.size(); i++)
	{
		Edge* v = model->m_edges[i];
		if (v->isScreenPointClose(pt, transform))
		{
			if (v->getDepth() < minDepth)
			{
				minDepth = v->getDepth();
				idSelected = i;
			}
		}
	}
	if (!(idSelected == -1))
	{
		Edge* selectedEdge = model->m_edges[idSelected];
		HalfEdge* validHed = m_faceBeingCollected->getValidHed(selectedEdge);
		if (!(validHed == nullptr))
		{
			selectedEdge->setColorEdge(Edge::selectedColor);
			m_faceBeingCollected->insertHed(validHed);
			return true;
		}
	}
	return false;
}

void Collector::insertNewVertex(Vertex* newVertex)
{
	m_hedBeingCollected->initPoint(newVertex);
}

bool Collector::isHedComplete()
{
	return m_hedBeingCollected->isComplete();
}

void Collector::reset()
{
	delete m_hedBeingCollected;
	delete m_faceBeingCollected;
	m_hedBeingCollected = nullptr;
	m_faceBeingCollected = nullptr;
	m_startCollecting = false;
}

bool Collector::isFaceComplete()
{
	return m_faceBeingCollected->isComplete();
}
