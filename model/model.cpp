#include <QtMath>
#include "model.h"

Model::Model()
{
	m_centroide = QVector3D(0.0, 0.0, 0.0);
	//// points
	//QVector3D* p0;
	//QVector3D* p1;
	//QVector3D* p2;
	//QVector3D* p3;

	//p3 = new QVector3D(0.0, 0.0, 1.0);
	//p1 = new QVector3D(1.0, 0.0, 0.0);
	//p2 = new QVector3D(qCos(2 * M_PI / 3), qSin(2 * M_PI / 3), 0.0);
	//p0 = new QVector3D(qCos(-2 * M_PI / 3), qSin(-2 * M_PI / 3), 0.0);





	//m_vertexes.push_back(new Vertex(p0,0));
	//m_vertexes.push_back(new Vertex(p1,1));
	//m_vertexes.push_back(new Vertex(p2,2));
	//m_vertexes.push_back(new Vertex(p3,3));



	//m_heds.push_back(new HalfEdge(m_vertexes[0], m_vertexes[1], 0));
	//m_heds.push_back(new HalfEdge(m_vertexes[1], m_vertexes[2], 1));
	//m_heds.push_back(new HalfEdge(m_vertexes[2], m_vertexes[0], 2));
	//m_heds.push_back(new HalfEdge(m_vertexes[0], m_vertexes[3], 3));
	//m_heds.push_back(new HalfEdge(m_vertexes[1], m_vertexes[3], 4));
	//m_heds.push_back(new HalfEdge(m_vertexes[2], m_vertexes[3], 5));

	//m_edges.push_back(new Edge(m_heds[0], 0));
	//m_edges.push_back(new Edge(m_heds[1], 1));
	//m_edges.push_back(new Edge(m_heds[2], 2));
	//m_edges.push_back(new Edge(m_heds[3], 3));
	//m_edges.push_back(new Edge(m_heds[4], 4));
	//m_edges.push_back(new Edge(m_heds[5], 5));


	////tri 0

	//m_faces.push_back(new Face (m_edges[0], m_edges[1], m_edges[2] , 0) );

	////tri1

	//m_faces.push_back(new Face(m_edges[0], m_edges[4], m_edges[3] , 1) );

	////tri2

	//m_faces.push_back(new Face(m_edges[1], m_edges[5], m_edges[4] , 2) );

	////tri3

	//m_faces.push_back(new Face(m_edges[2], m_edges[3], m_edges[5], 3) );












}

void Model::insertNewHed(HalfEdge* hed)
{
	int nHed = m_heds.size();
	// create new hed 
	HalfEdge* newHed = new HalfEdge(hed->getP1(), hed->getP2(), nHed);
	HalfEdge* hedTwin = new HalfEdge(hed->getP2(), hed->getP1(), nHed + 1);
	m_heds.push_back(newHed);
	m_heds.push_back(hedTwin);
	// create also a new Edge
	Edge* newEdge = new Edge(newHed, hedTwin, m_edges.size());
	m_edges.push_back(newEdge);
	// update heds
	newHed->setEdge(newEdge);
	hedTwin->setEdge(newEdge);
	


}

void Model::insertNewVertex(QVector3D* pt)
{
	int count = m_vertexes.size();
	m_vertexes.push_back( new Vertex(*pt , count) );
	updateCentroide();
}

void Model::insertNewFace(Face* newFace)
{
	int n = m_faces.size();
	std::vector<HalfEdge*> heds = newFace->getHeds();
	int nHeds = heds.size();
	switch (nHeds)
	{
	case 3:
	{
		// one triangle
		Face* newFace = new Face(heds[0], heds[1], heds[2], n);
		m_faces.push_back(newFace);
		break;
	}
	case 4:
	{
		// two triangles
		// first create a new hed dividing the triangles
		HalfEdge* newHed = new HalfEdge(heds[0]->getP1(), heds[1]->getP2(), m_heds.size());
		HalfEdge* newHedTwin = new HalfEdge(heds[1]->getP2(), heds[0]->getP1(), m_heds.size() + 1);
		Edge* newEdge = new Edge(newHed, newHedTwin, m_edges.size());
		newHed->setEdge(newEdge);
		newHedTwin->setEdge(newEdge);
		m_heds.push_back(newHed);
		m_heds.push_back(newHedTwin);
		m_edges.push_back(newEdge);
		// build two triangles
		Face* newFace1 = new Face(heds[0], heds[1], newHedTwin, n);
		Face* newFace2 = new Face(newHed, heds[2], heds[3], n + 1);
		m_faces.push_back(newFace1);
		m_faces.push_back(newFace2);

		break;
	}
	}
}

QVector3D* Model::getCentroide()
{
	return &m_centroide;
}

void Model::updateCentroide()
{
	int n = m_vertexes.size();
	Vertex* lastVertex = m_vertexes.back();
	m_centroide = ((n - 1) * m_centroide + *lastVertex->getP() ) / n;
}

void Model::updateIds()
{
	for (int i = 0; i < m_faces.size(); i++)
	{
		 m_faces[i]->setId(i);
	}
	for (int i = 0; i < m_vertexes.size(); i++)
	{
		m_vertexes[i]->setId(i);
	}
	for (int i = 0; i < m_edges.size(); i++)
	{
		m_edges[i]->setId(i);
	}

}

void Model::deleteSO()
{
	
	bool faceRemoved , vertexRemoved, edgeRemoved;
	do {
		faceRemoved = false;
		for (int i = 0; i < m_faces.size(); ++i) {
			Face* face = m_faces[i];
			if (face->isSelected()) {
				deleteFace(i);
				faceRemoved = true;
				break;
			}
		}
	} while (faceRemoved);

	do {
		edgeRemoved = false;
		for (int i = 0; i < m_edges.size(); ++i) {
			Edge* edge = m_edges[i];
			if (edge->isSelected()) {
				// delete also other objects binded
				do {
					faceRemoved = false;
					for (int j = 0; j < m_faces.size(); ++j) {
						Face* face = m_faces[j];
						if (face->contain(edge)) {
							deleteFace(j);
							faceRemoved = true;
							break;
						}
					}
				} while (faceRemoved);
				deleteEdge(i);
				edgeRemoved = true;
				break;
			}
		}
	} while (edgeRemoved);

	do {
		vertexRemoved = false;
		for (int i = 0; i < m_vertexes.size(); ++i) {
			Vertex* v = m_vertexes[i];
			if (v->isSelected()) {
				// delete also the other objects binded

				do {
					edgeRemoved = false;
					for (int j = 0; j < m_edges.size(); ++j) {
						Edge* edge = m_edges[j];
						if (edge->contain(v)) {
							// delete also other objects binded
							do {
								faceRemoved = false;
								for (int k = 0; k < m_faces.size(); ++k) {
									Face* face = m_faces[k];
									if (face->contain(edge)) {
										deleteFace(k);
										faceRemoved = true;
										break;
									}
								}
							} while (faceRemoved);
							deleteEdge(j);
							edgeRemoved = true;
							break;
						}
					}
				} while (edgeRemoved);

				deleteVertex(i);
				vertexRemoved = true;
				break;
			}
		}
	} while (vertexRemoved);

}

void Model::deleteVertex(int pos)
{
	delete m_vertexes[pos];
	m_vertexes.erase(m_vertexes.begin() + pos);
}

void Model::deleteEdge(int pos)
{
	Edge* edge =  m_edges[pos];
	HalfEdge* hed = edge->getCurrentHed();
	HalfEdge* twin = edge->getTwinHed();
	m_edges.erase(m_edges.begin() + pos);

	m_heds.erase(m_heds.begin() + hed->getId());
	m_heds.erase(m_heds.begin() + hed->getId());
	delete edge, hed, twin;
	// update heds
	for (int i = 0; i < m_heds.size(); i++)
	{
		m_heds[i]->setId(i);
		
	}
}

void Model::deleteFace(int pos)
{
	delete m_faces[pos];
	m_faces.erase(m_faces.begin() + pos);
}


