#ifndef VERTEX_H
#define VERTEX_H

#include <QVector3D>
#include <QVector2D>
#include <QVector4D>
#include <vector>

namespace topoApp {
	class Vertex
	{
	public:
		Vertex();
		Vertex(QVector3D& p, int id);
		QVector3D* getP();
		QVector3D getColor();
		QVector3D getSelectedColor();
		QVector3D getNormalColor();
		qreal getDepth();
		std::vector<QVector3D>* getBB();
		void buildBB();
		void setColor(QVector3D color);
		void getBBOnClippingCoordinates(qreal* xmin, qreal* xmax, qreal* ymin, qreal* ymax, QMatrix4x4* transform);
		bool isScreenPointClose(const QVector2D point, QMatrix4x4* transform);
		int getId();
		qreal X();
		qreal Y();
		qreal Z();
		void setPosValue(int pos, qreal value);
		// fixed color 
		static const QVector3D normalColor;
		static const QVector3D selectedColor;
		void setIsSelected(bool select) { m_isSelected = select; };
		bool isSelected() { return m_isSelected; };
		void setId(int id) { m_id = id; };
		bool getBC() { return m_isBCActive; };
		void setBC(bool bc) { m_isBCActive = bc; };
		qreal getTemperature() { return m_temperature; };
		void setTemperature(qreal temperature) { m_temperature = temperature; };
	private:
		//bounding box
		QVector3D m_p;
		std::vector<QVector3D> m_bb;
		QVector3D m_color; // current color
		const qreal m_radius = 0.055;

		qreal m_depth;
		int m_id;
		bool m_isSelected = false;
		bool m_isBCActive = false;
		qreal m_temperature;


	};


}
#endif // !POINT_H
