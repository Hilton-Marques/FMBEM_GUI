#ifndef GLCANVAS_H
#define GLCANVAS_H
#include <QQuaternion>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QSurfaceFormat>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLShader>

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QMatrix4x4>
#include <QMatrix2x2>
#include <QQuaternion>
#include <QVector2D>
#include <QBasicTimer>
#include <QOpenGLShaderProgram>
#include <QOpenGLTexture>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include <QString>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QtWidgets/QTableWidget>
#include <QTableWidgetItem>

#include "model/model.h"
#include "surfaces/sphere.h"
#include "surfaces/cylinder.h"
#include "surfaces/cone.h"
#include "surfaces/arrow.h"

#include "canvas/handleselectedobjects.h"
#include "canvas/collector.h"

using namespace topoApp;

enum class ActionType {
	UNDEF, SELECTION, COLLECTION
};

class GLcanvas : public QOpenGLWidget, protected QOpenGLFunctions
{
	Q_OBJECT
public:
	GLcanvas(QWidget* parent = 0);
	void setModel(Model* model);
	void setLogTables(QTableWidget* logPointTable, QTableWidget* logEdgesTable, QTableWidget* logFacesTable);
	void scaleWordWindow(qreal fac);
	void fit2world();
	void setMouseAction(ActionType actionType);
	Collector m_collector;
	HandleSO* m_handleSO;
	void setNewEyeReference(QVector3D* ref);
	void deleteSelectedObjects();
	void reset();
protected:
	void initializeGL();
	void paintGL();
	void resizeGL(int w, int h);
	void initShaders(QOpenGLShaderProgram* shader, QString name);
	void initModel();
	void initSphere();
	void initCylinder();
	void initCone();
	void initArrow();

	void panEyeWindow(QVector2D panV);
	
	void mousePressEvent(QMouseEvent* event) override;
	void mouseMoveEvent(QMouseEvent* event) override;
	void mouseReleaseEvent(QMouseEvent* event) override;
	void wheelEvent(QWheelEvent* event);
	void keyPressEvent(QKeyEvent* event) override;
	void keyReleaseEvent(QKeyEvent* event) override;
	QVector2D convertRasterPtToEye(QVector2D rasterPt);
	void rotateViewAboutX(qreal mouseMoveFacY);
	void rotateViewAboutY(qreal mouseMoveFacX);
	QMatrix4x4 rotateAboutAfixedPoint(qreal angle, QVector3D axis, QVector3D point);
	void buildEyeSystem();
	QMatrix4x4 m_projectionMatrix;
	void checkFaceOrientation(Face* face);
	
	

private slots:
	void on_logPoints_itemChanged(QTableWidgetItem* item);
	void on_logPoints_doubleClicked(QTableWidgetItem* item);

	void on_logEdges_itemChanged(QTableWidgetItem* item);
	void on_logEdges_doubleClicked(QTableWidgetItem* item);

	void on_logFaces_itemChanged(QTableWidgetItem* item);
	void on_logFaces_doubleClicked(QTableWidgetItem* item);




private:

	Model* m_model;
	QOpenGLShaderProgram* m_shaderProgram;
	int m_mouseMoveTol = 2;

	QOpenGLVertexArrayObject* vao;
	QOpenGLVertexArrayObject* sphereVao;
	QOpenGLVertexArrayObject* cylinderVao;
	QOpenGLVertexArrayObject* coneVao;
	QOpenGLVertexArrayObject* arrowVao;


	QOpenGLBuffer m_modelArray;
	QOpenGLBuffer m_modelIndex;

	QOpenGLBuffer m_sphereArray;
	QOpenGLBuffer m_sphereIndex;

	QOpenGLBuffer m_cylinderArray;
	QOpenGLBuffer m_cylinderIndex;

	QOpenGLBuffer m_coneArray;
	QOpenGLBuffer m_coneIndex;

	QOpenGLBuffer m_arrowArray;
	QOpenGLBuffer m_arrowIndex;

	// camera position 
	QVector3D m_eye ;
	QVector3D m_ref;
	QVector3D m_up;
	// Matrixes
	QMatrix4x4 m_modelMatrix;
	QMatrix4x4 m_viewMatrix;

	QMatrix4x4 m_mvp;

	Sphere m_sphere;
	Cylinder m_cylinder;
	Cone m_cone;
	Arrow m_arrow;

	// mouse events
	ActionType m_currentAction = ActionType::SELECTION;
	bool m_buttonPressed = false;
	//bool m_isSelecting = true;
	QVector2D m_pt0;
	QVector2D m_pt1;
	Qt::MouseButton m_mouseButton;
	// bounding box
	double m_left, m_right, m_bottom, m_top, m_front, m_back;
	double m_leftS, m_rightS, m_bottomS, m_topS, m_frontS, m_backS; // static bounding box (pan and zoom could change these coordinates)
	//camera system
	QVector3D m_xEye, m_yEye, m_zEye;
	GLint m_mvpLoc, m_colorLoc, m_mv, m_cofMv;
	// hidden faces
	std::vector<Face*> m_hiddenFaces;
	std::vector<Face*> m_staticFaces;
	// Log tables
	QTableWidget* m_logPoints;
	QTableWidget* m_logEdges;
	QTableWidget* m_logFaces;

	// max size of screen
	const qreal m_max_size = 10.0;
	
};



#endif // !GLCANVAS_H
