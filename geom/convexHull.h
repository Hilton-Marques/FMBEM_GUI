#ifndef CONVEX_HULL_H
#define CONVEX_HULL_H
#include <QtMath>
#include <QVector3D>
#include <QMatrix4x4>

#include "../model/model.h"
#include "../topoApp/VertexApp.h"
#include "../topoApp/FaceApp.h"
#include "../topoApp/HalfEdgeApp.h"
#include "../topoApp/EdgeApp.h"

using namespace topoApp;

#include <vector>
class ConvexHull
{
public:
	ConvexHull(Model* model);
	void initModel();
	Face findTriangleOnHull();
	Vertex* pivotOnEdge(HalfEdge* hed);
	HalfEdge* findEdgeOnHull();
	static qreal squaredArea(QVector3D* p0, QVector3D *p1, QVector3D *p2);
	static qreal signedVolume(QVector3D* p0, QVector3D* p1, QVector3D* p2, QVector3D* p3);

private:
	std::vector<Vertex*> m_verts;
	Model* m_model;

};






#endif // !CONVEX_HULL_H
