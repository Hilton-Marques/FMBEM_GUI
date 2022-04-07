#ifndef EDGE_H
#define EDGE_H

#include <QVector3D>
#include <QVector2D>

#include "HalfEdgeApp.h"
#include "VertexApp.h"

namespace topoApp {

	class Edge
	{
	public:
		Edge();
		Edge(HalfEdge* hed, HalfEdge* twin, int id);

		QVector3D* p1();
		QVector3D* p2();

		void buildBB();
		void getBBOnClippingCoordinates(qreal* xmin, qreal* xmax, qreal* ymin, qreal* ymax, QMatrix4x4* transform);
		qreal getDepth();

		bool isScreenPointClose(QVector2D point, QMatrix4x4* transform);

		QVector3D getColorEdge();
		QVector3D getColorCone();
		QVector3D getNormalColorEdge();
		QVector3D getNormalColorCone();

		QVector3D getSelectedColorEdge();
		QVector3D getSelectedColorCone();
		void setColorEdge(QVector3D color);
		void setColorCone(QVector3D color);
		void createLocalCoordinateSystem();
		int getId();
		int getIdP1();
		int getIdP2();
		void setNewPoint(int pos, Vertex* point);
		HalfEdge* getCurrentHed();
		HalfEdge* getTwinHed();
		HalfEdge* getTwinHed(HalfEdge* hed);

		static const QVector3D normalColor;
		static const QVector3D selectedColor;
		void setIsSelected(bool select) { m_isSelected = select; };
		bool isSelected() { return m_isSelected; };
		bool contain(Vertex* v);
		void setId(int id) { m_id = id; };
		void setIsProcessed(bool processed) { m_isProcessed = processed; };
		bool isProcessed() { return m_isProcessed; };
	private:
		qreal m_radius = 0.012;

		HalfEdge* m_hed;
		HalfEdge* m_twinHed;

		Vertex* m_p1;
		Vertex* m_p2;
		std::vector<QVector3D> m_bb;
		QVector3D m_centroide;
		qreal m_depth;
		QVector3D m_colorEdge; // current color
		QVector3D m_colorCone; // current color
		//create a local coordinate system
		QVector3D m_x, m_y, m_z;
		int m_id;
		bool m_isSelected = false;
		bool m_isProcessed = false;
	};
}
#endif
