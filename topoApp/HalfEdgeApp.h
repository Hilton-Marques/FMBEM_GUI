#ifndef HALFEDGE_H
#define HALFEDGE_H

#include "VertexApp.h"

namespace topoApp {
	class Edge;

	class HalfEdge
	{
	public:
		HalfEdge();
		HalfEdge(Vertex* p1, Vertex* p2, int id);
		Vertex* getP1();
		Vertex* getP2();
		int getId();
		void setNewPoint(int pos, Vertex* value);
		void initPoint(Vertex* newpt);
		bool isComplete();
		void setId(int id);
		void setEdge(Edge* edge) { m_edge = edge; };
		Edge* getEdge() { return m_edge; };
		HalfEdge* getTwin();
	private:
		Vertex* m_p1 = nullptr;
		Vertex* m_p2 = nullptr;
		int m_id;
		Edge* m_edge ;
	};

}
#endif
