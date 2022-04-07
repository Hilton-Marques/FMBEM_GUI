#ifndef HANDLESELECTEDOBJECTS_H
#define HANDLESELECTEDOBJECTS_H
#include <vector>
#include <QMatrix4x4>
#include <QVector2D>
#include <QString>
#include <QStringList>
#include <QtWidgets/QTableWidget>
#include <QTableWidgetItem>

#include "../model/model.h"
#include "../topoApp/VertexApp.h"
#include "../topoApp/EdgeApp.h"
#include "../topoApp/FaceApp.h"

using namespace topoApp;
enum class SelectionType {
	SELECTION, SELECTIONBC
};

class HandleSO
{
public:
	HandleSO();
	HandleSO(Model* model, QTableWidget* logPointTable,QTableWidget* logEdgeTable, QTableWidget* logFaceTable);

	void selectObjects(QVector2D point, QMatrix4x4* transform);
	void selectVertices(QVector2D point, QMatrix4x4* transform);
	void selectEdges(QVector2D point, QMatrix4x4* transform);
	void selectFaces(QVector2D point, QMatrix4x4* transform);
	std::vector<Face*> getSelectedFaces();
	void reset();
	void changeColorSO();
	bool isAnySelected();
	bool m_isPointSelected = false;
	bool m_isEdgeSelected = false;
	bool m_isFaceSelected = false;

	bool m_logPointsInit = false;
	bool m_logEdgesInit = false;
	bool m_logFacesInit = false;

	void updateTable();
	void deleteSelectedObjects();
	void setSelectionType(SelectionType selection) { m_selectionType = selection; };
	SelectionType getSelectionType() { return m_selectionType;  };
private:
	std::vector<Vertex*> m_selectedPoints;
	std::vector<Edge*> m_selectedEdges;
	std::vector<Face*> m_selectedFaces;
	Model* m_model;
	QTableWidget* m_logPoints;
	QTableWidget* m_logEdges;
	QTableWidget* m_logFaces;
	SelectionType m_selectionType = SelectionType::SELECTION;
	// handle talble inits
	


};

#endif // !HANDLESELECTEDOBJECTS_H
