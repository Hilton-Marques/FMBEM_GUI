#ifndef MODEL_H
#define MODEL_H
#include <QOpenGLFunctions>
#include <vector>
#include <QVector3D>
#include "../topoApp/VertexApp.h"
#include "../topoApp/EdgeApp.h"
#include "../topoApp/HalfEdgeApp.h"
#include "../topoApp/FaceApp.h"

using namespace topoApp;
class Model
{
public:
	Model();
	std::vector<Vertex*> m_vertexes;
	std::vector<Edge*> m_edges;
	std::vector<Face*> m_faces;
	std::vector<HalfEdge*> m_heds;
	void insertNewHed(HalfEdge* newHed);
	void insertNewVertex(QVector3D* newVertex);
	void insertNewFace(Face* newFace);
	QVector3D* getCentroide();
	void updateCentroide();
	void updateIds();
	void deleteSO();
	void deleteVertex(int pos);
	void deleteEdge(int pos);
	void deleteFace(int pos);
private:
	QVector3D m_centroide;


};
#endif // !MODEL_H
