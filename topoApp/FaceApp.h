#ifndef FACE_H
#define FACE_H
#include <QVector4D>
#include <QVector3D>
#include <QVector2D>
#include "EdgeApp.h"
#include <vector>

namespace topoApp {
	class Face
	{
	public:
		Face();
		Face(HalfEdge* hed1, HalfEdge* hed2, HalfEdge* hed3, int id);
		QVector3D getCentroide();
		std::vector<QVector3D> getPoints();
		bool isScreenPointClose(QVector2D point, QMatrix4x4* transform);
		QVector3D getNormalColor();
		QVector3D getSelectedColor();
		QVector3D getCurrentColor();
		void setColor(QVector3D color);
		qreal getDepth();
		static qreal orientedArea(QVector2D a, QVector2D b, QVector2D c);
		int getId();
		void setId(int id);
		Edge* getE1();
		Edge* getE2();
		Edge* getE3();
		HalfEdge* getValidHed(Edge* edge);
		void insertHed(HalfEdge* hed);
		bool isComplete();
		void reset();
		std::vector<HalfEdge*> getHeds() { return m_heds; }; // could be more than 3
		void setIsSelected(bool select) { m_selected = select; };
		bool isSelected() { return m_selected; };
		bool contain(Edge* edge);
		bool getBC() { return m_isBCActive; };
		void setBC(const bool bc) { m_isBCActive = bc; };
		qreal getGradient() { return m_gradient; };
		void setGradient(qreal gradient) { m_gradient = gradient; };
		QVector3D getNormal();
		void reverseHeds();
		qreal getArea();
		std::vector<HalfEdge*> get3Heds();
		std::vector<QVector3D> getPointsAndNormals();
	private:
		HalfEdge* m_he1;
		HalfEdge* m_he2;
		HalfEdge* m_he3;
		QVector3D* m_p1;
		QVector3D* m_p2;
		QVector3D* m_p3;
		QVector3D m_color;
		qreal m_depth;
		bool m_isHidden = false;
		int m_id;
		std::vector<HalfEdge*> m_heds;
		bool m_selected = false;
		bool m_isBCActive = true;
		qreal m_gradient = 1.0;
	};

}
#endif // !FACE_H
