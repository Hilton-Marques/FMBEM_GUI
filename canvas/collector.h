#ifndef COLLECTOR_H
#define COLLECTOR_H

#include <QVector2D>
#include <QMatrix4x4>
#include "../topoApp/HalfEdgeApp.h"
#include "../topoApp/FaceApp.h"
#include "../model/model.h"

using namespace topoApp;

enum class ObjectType {
	UNDEF, HALFEDGE, FACE
};

class Collector
{
public:
	Collector();
	~Collector();
  void startCollection();
	void setObjectType(ObjectType type);
	ObjectType getObjectType();
	bool isPointSelected(QVector2D& pt, QMatrix4x4* transform, const Model* model);
	bool isEdgeSelected(QVector2D &pt, QMatrix4x4* transform,const Model* model);
	void insertNewVertex(Vertex* newVertex);
	bool isHedComplete();
	HalfEdge* getCollectedHed() { return m_hedBeingCollected; };
	void reset();
	bool m_startCollecting = false;
	Face* getCollectedFace() { return m_faceBeingCollected; };
	bool isFaceComplete();
private:
	ObjectType m_objectType = ObjectType::UNDEF;
	HalfEdge* m_hedBeingCollected = nullptr;
	Face* m_faceBeingCollected = nullptr;
};


#endif