#include "glcanvas.h"

#include <QGuiApplication>
#include <QMatrix4x4>
#include <QOpenGLShaderProgram>
#include <QScreen>
#include <QtMath>

#include <QMouseEvent>
#include <QKeyEvent>
#include <algorithm>    // std::remove


GLcanvas::GLcanvas(QWidget* parent):QOpenGLWidget(parent), m_modelIndex(QOpenGLBuffer::IndexBuffer)
{

	QSurfaceFormat format = QSurfaceFormat::defaultFormat();
	format.setSamples(16);
	format.setSwapInterval(0);
	QSurfaceFormat::setDefaultFormat(format);
	this->setFormat(format);
	// camera
	m_eye = QVector3D(0.0, 0.0, 8.0);
	m_ref = QVector3D(0.0, 0.0, 0.0);
	m_up = QVector3D(0.0, 1.0, 0.0);
	// bounding box
	m_left = -2.0;
	m_right = 2.0;
	m_bottom = -2.0;
	m_top = 2.0;
	m_front = 5.0;
	m_back = 15.0;

	

	
}

void GLcanvas::initializeGL()
{

	initializeOpenGLFunctions();
	setFocusPolicy(Qt::StrongFocus);

	m_shaderProgram = new QOpenGLShaderProgram;



	initShaders(m_shaderProgram, "face");

	m_mvpLoc = m_shaderProgram->uniformLocation("mvp");
	m_colorLoc = m_shaderProgram->uniformLocation("color");
	m_mv = m_shaderProgram->uniformLocation("mv");
	m_cofMv = m_shaderProgram->uniformLocation("cofMv");

	initModel();
	initSphere();
	initCylinder();
	initCone();
	initArrow();
	

}

void GLcanvas::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
	//glCullFace(GL_BACK);//GL_BACK
	glClearColor(0.2,0.2, 0.2, 1.0);

	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//Initialize Matrixes
	
	QMatrix4x4 matrix, model, camera,rot, mv, cofMv;
	camera.lookAt(m_eye, m_ref, m_up);

	//to Cylinder and Cone transformations
	QVector3D u, v, perp;
	u = QVector3D(0.0, 0.0, 1.0);

	m_shaderProgram->bind();

	//
	vao->bind();
	m_modelArray.bind();

	model.setToIdentity();
	mv.setToIdentity();
	cofMv.setToIdentity();
	matrix = m_projectionMatrix * camera * model;
	mv = camera * model;
	cofMv = (mv.inverted()).transposed();

	m_shaderProgram->setUniformValue(m_mvpLoc, matrix);
	m_shaderProgram->setUniformValue(m_mv, mv);
	m_shaderProgram->setUniformValue(m_cofMv, cofMv);

	int count = 0;
	for (Face* face : m_model->m_faces)
	{
		m_modelArray.allocate(face->getPointsAndNormals().data(), 6 * sizeof(QVector3D));
		m_shaderProgram->setUniformValue(m_colorLoc, face->getCurrentColor());
		glDrawArrays(GL_TRIANGLES, 0, 3);
	}
	//glDrawElements(GL_TRIANGLES, m_model->m_indices.size(), GL_UNSIGNED_SHORT, 0);

	m_modelArray.release();
	m_modelIndex.release();
	vao->release();

	sphereVao->bind();
	m_sphereArray.bind();
	m_sphereIndex.bind();

	for (Vertex* v : m_model->m_vertexes)
	{
		model.setToIdentity();
		mv.setToIdentity();
		cofMv.setToIdentity();

		model.translate(*v->getP());
		matrix = m_projectionMatrix * camera * model;
		
		mv = camera * model;
		cofMv = (mv.inverted()).transposed();

		m_shaderProgram->setUniformValue(m_mvpLoc, matrix);
		m_shaderProgram->setUniformValue(m_mv, mv);
		m_shaderProgram->setUniformValue(m_cofMv, cofMv);

		m_shaderProgram->setUniformValue(m_colorLoc, v->getColor());
		glDrawElements(GL_QUADS, m_sphere.m_indices.size(), GL_UNSIGNED_SHORT, 0);
	}
	if (m_handleSO->getSelectionType() == SelectionType::SELECTIONBC)
	{
		for (Vertex* v : m_model->m_vertexes)
		{
			if (v->getBC())
			{
				model.setToIdentity();
				mv.setToIdentity();
				cofMv.setToIdentity();

				model.translate(*v->getP());
				model.scale(1.1, 1.1, 1.1);
				matrix = m_projectionMatrix * camera * model;

				mv = camera * model;
				cofMv = (mv.inverted()).transposed();

				m_shaderProgram->setUniformValue(m_mvpLoc, matrix);
				m_shaderProgram->setUniformValue(m_mv, mv);
				m_shaderProgram->setUniformValue(m_cofMv, cofMv);

				m_shaderProgram->setUniformValue(m_colorLoc, QVector3D(1.0,0.8,0.0));
				glDrawElements(GL_QUADS, m_sphere.m_indices.size(), GL_UNSIGNED_SHORT, 0);

			}

		}

	}


	m_sphereArray.release();
	m_sphereIndex.release();
	sphereVao->release();


	cylinderVao->bind();
	m_cylinderArray.bind();
	m_cylinderIndex.bind();
	

	
	for (Edge* e : m_model->m_edges )
	{
		v = *e->p2() - *e->p1();
		qreal d = v.length();
		perp = QVector3D::crossProduct(u, v.normalized());
		qreal angle = qAcos(QVector3D::dotProduct(u, v.normalized()));
		angle = qRadiansToDegrees(angle);

		model.setToIdentity();
		mv.setToIdentity();
		cofMv.setToIdentity();

		model.translate(0.5 * (*e->p2() + *e->p1()));
		model.rotate(angle, perp);
		model.scale(1.0, 1.0, d);
		
		matrix = m_projectionMatrix * camera * model;


		mv = camera * model;
		cofMv = (mv.inverted()).transposed();

		m_shaderProgram->setUniformValue(m_mvpLoc, matrix);
		m_shaderProgram->setUniformValue(m_mv, mv);
		m_shaderProgram->setUniformValue(m_cofMv, cofMv);


		m_shaderProgram->setUniformValue(m_colorLoc, e->getColorEdge());

		glDrawElements(GL_QUADS, m_cylinder.m_indices.size(), GL_UNSIGNED_SHORT, 0);

	}

	m_cylinderArray.release();
	m_cylinderIndex.release();

	cylinderVao->release();

	//Axis
	QVector3D x(1.0, 0.0, 0.0);
	QVector3D y(0.0, 1.0, 0.0);
	QVector3D z(0.0, 0.0, 1.0);
	std::vector<QVector3D> axis = { x,y,z };

	////Cone

	//coneVao->bind();
	//m_coneArray.bind();
	//m_coneIndex.bind();



	//m_shaderProgram->setUniformValue(m_colorLoc, QVector3D(0.6,1.0,0.0));

	//for (Edge* e : m_model->m_edges)
	//{

	//	v = *(e->p2()) - *(e->p1());
	//	qreal d = v.length();
	//	perp = QVector3D::crossProduct(u, v.normalized());
	//	qreal angle = qAcos(QVector3D::dotProduct(u, v.normalized()));
	//	angle = qRadiansToDegrees(angle);

	//	model.setToIdentity();
	//	model.translate(0.5 * ( *(e->p2()) + *(e->p1()) ) );
	//	model.rotate(angle, perp);
	//	model.scale(0.05, 0.05, 0.2);


	//	matrix = m_projectionMatrix * camera * model;
	//	m_shaderProgram->setUniformValue(m_mvpLoc, matrix);
	//	m_shaderProgram->setUniformValue(m_colorLoc, e->getColorCone());

	//	glDrawElements(GL_QUADS, m_cone.m_indices.size(), GL_UNSIGNED_SHORT, 0);

	//}

	//m_coneArray.release();
	//m_coneIndex.release();
	//coneVao->release();


	arrowVao->bind();
	m_arrowArray.bind();
	m_arrowIndex.bind();

	for (QVector3D axi : axis)
	{
		perp = QVector3D::crossProduct(u, axi);
		qreal angle = qAcos(QVector3D::dotProduct(u, axi));
		angle = qRadiansToDegrees(angle);
		model.setToIdentity();
		mv.setToIdentity();
		cofMv.setToIdentity();
		model.rotate(angle, perp);

		matrix = m_projectionMatrix * camera * model;
		mv = camera * model;
		cofMv = (mv.inverted()).transposed();

		m_shaderProgram->setUniformValue(m_mvpLoc, matrix);
		m_shaderProgram->setUniformValue(m_mv, mv);
		m_shaderProgram->setUniformValue(m_cofMv, cofMv);

		m_shaderProgram->setUniformValue(m_colorLoc, axi);

		glDrawElements(GL_QUADS, m_arrow.m_indices.size(), GL_UNSIGNED_SHORT, 0);


	}
	if (m_handleSO->getSelectionType() == SelectionType::SELECTIONBC)
	{
		QVector3D normal, centroide, translation;
		QVector3D color(1.0, 0.8, 0.0);
		QMatrix4x4 rotate, translate,scale;
		qreal area;
		for (Face* face : m_model->m_faces)
		{
			if (face->getBC())
			{
				normal = face->getNormal();
				centroide = face->getCentroide();
				area = face->getArea();
				perp = QVector3D::crossProduct(u, normal);
				qreal angle = qAcos(QVector3D::dotProduct(u, normal));
				angle = qRadiansToDegrees(angle);
				rotate.setToIdentity();
				rotate.rotate(angle, perp);
				scale.setToIdentity();
				qreal value = face->getGradient() * 0.5;
				scale.scale(1.0, 1.0, value);

				QVector3D p1, p2, mid, x,y;
				qreal radiusTemp ;
				qreal radius = FLT_MAX;
				for (HalfEdge* hed : face->get3Heds())
				{
					p1 = *hed->getP1()->getP();
					p2 = *hed->getP2()->getP();
					mid = 0.5 * (p1 + p2);
					radiusTemp = (mid - centroide).length();
					if (radiusTemp < radius)
					{
						radius = radiusTemp;
						x = (mid - centroide).normalized();
					}
				}
				y = QVector3D::crossProduct(normal, x);
				int s = 20;
				qreal fac = 0.7;
				qreal h = 2*M_PI / s;
				qreal cosTeta, sinTeta;
				for (int i = 0; i < s; i++)
					{
					cosTeta = std::cos(i * h);
					sinTeta = std::sin(i * h);
					translate.setToIdentity();
					translation = centroide + fac*radius*cosTeta*x + fac*radius*sinTeta*y;
					translate.translate(translation);
					model.setToIdentity();
					mv.setToIdentity();
					cofMv.setToIdentity();

					model = translate * rotate * scale;
					mv = camera * model;
					cofMv = (mv.inverted()).transposed();

					matrix = m_projectionMatrix * camera * model;

					m_shaderProgram->setUniformValue(m_mvpLoc, matrix);
					m_shaderProgram->setUniformValue(m_mv, mv);
					m_shaderProgram->setUniformValue(m_cofMv, cofMv);
					m_shaderProgram->setUniformValue(m_colorLoc, color);
					glDrawElements(GL_QUADS, m_arrow.m_indices.size(), GL_UNSIGNED_SHORT, 0);
				}
				//translate.setToIdentity();
				//translate.translate(centroide);
				//model.setToIdentity();
				//model = translate * rotate * scale;
				//matrix = m_projectionMatrix * camera * model;
				//m_shaderProgram->setUniformValue(m_mvpLoc, matrix);
				//m_shaderProgram->setUniformValue(m_colorLoc, color);
				//glDrawElements(GL_QUADS, m_arrow.m_indices.size(), GL_UNSIGNED_SHORT, 0);
			
			}

		}

	}

	m_arrowArray.release();
	m_arrowIndex.release();
	arrowVao->release();



	m_shaderProgram->release();

}

void GLcanvas::resizeGL(int w, int h)
{
	glViewport(0, 0, w, h);
	if (!m_buttonPressed) //when a click a point, the table changes so the resize is called
	{
		fit2world();
	}
	
}

void GLcanvas::initShaders(QOpenGLShaderProgram* shader, QString name)
{
	QString filenameVert = name + ".vert";
	QString filenameFrag = name + ".frag";

	if (!shader->addShaderFromSourceFile(QOpenGLShader::Vertex, filenameVert))
		close();

	if (!shader->addShaderFromSourceFile(QOpenGLShader::Fragment, filenameFrag))
		close();

	if (!shader->link())
		close();

	if (!shader->bind())
		close();
}

void GLcanvas::initModel()
{
	vao = new QOpenGLVertexArrayObject(this);
	vao->create();
	vao->bind();
	m_modelArray = QOpenGLBuffer(QOpenGLBuffer::VertexBuffer);
	m_modelArray.create();
	m_modelArray.setUsagePattern(QOpenGLBuffer::DynamicDraw);
	m_modelArray.bind();  
	m_shaderProgram->enableAttributeArray(0);
	m_shaderProgram->setAttributeBuffer(0, GL_FLOAT, 0, 3, 0);
	m_shaderProgram->enableAttributeArray(1);
	m_shaderProgram->setAttributeBuffer(1, GL_FLOAT, 3 * sizeof(QVector3D), 3, 0);
	m_modelArray.release();

	
	vao->release();
}

void GLcanvas::initSphere()
{
	sphereVao = new QOpenGLVertexArrayObject(this);
	sphereVao->create();
	sphereVao->bind();
	m_sphereArray = QOpenGLBuffer(QOpenGLBuffer::VertexBuffer);
	m_sphereArray.create();
	m_sphereArray.setUsagePattern(QOpenGLBuffer::StaticDraw);
	m_sphereArray.bind();
	m_sphereArray.allocate(m_sphere.m_vertices.data(), m_sphere.m_vertices.size() * sizeof(QVector3D));
	m_sphereIndex = QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);
	m_sphereIndex.create();
	m_sphereIndex.setUsagePattern(QOpenGLBuffer::StaticDraw);
	m_sphereIndex.bind();
	m_sphereIndex.allocate(m_sphere.m_indices.data(), m_sphere.m_indices.size() * sizeof(GLushort));
	m_shaderProgram->enableAttributeArray(0);
	m_shaderProgram->setAttributeBuffer(0, GL_FLOAT, 0, 3, 0);
	m_shaderProgram->enableAttributeArray(1);
	m_shaderProgram->setAttributeBuffer(1, GL_FLOAT, m_sphere.m_indices.size() * sizeof(QVector3D), 3, 0);
	m_sphereArray.release();
	m_sphereIndex.release();
	sphereVao->release();

}
void GLcanvas::initCylinder()
{
	cylinderVao = new QOpenGLVertexArrayObject(this);
	cylinderVao->create();
	cylinderVao->bind();
	m_cylinderArray = QOpenGLBuffer(QOpenGLBuffer::VertexBuffer);
	m_cylinderArray.create();
	m_cylinderArray.setUsagePattern(QOpenGLBuffer::StaticDraw);
	m_cylinderArray.bind();
	m_cylinderArray.allocate(m_cylinder.m_vertices.data(), m_cylinder.m_vertices.size() * sizeof(QVector3D));
	m_cylinderIndex = QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);
	m_cylinderIndex.create();
	m_cylinderIndex.setUsagePattern(QOpenGLBuffer::StaticDraw);
	m_cylinderIndex.bind();
	m_cylinderIndex.allocate(m_cylinder.m_indices.data(), m_cylinder.m_indices.size() * sizeof(GLushort));
	m_shaderProgram->enableAttributeArray(0);
	m_shaderProgram->setAttributeBuffer(0, GL_FLOAT, 0, 3, 0);
	m_shaderProgram->enableAttributeArray(1);
	m_shaderProgram->setAttributeBuffer(1, GL_FLOAT, m_cylinder.m_indices.size() * sizeof(QVector3D), 3, 0);
	m_cylinderArray.release();
	m_cylinderArray.release();
	cylinderVao->release();

}

void GLcanvas::initCone()
{
	coneVao = new QOpenGLVertexArrayObject(this);
	coneVao->create();
	coneVao->bind();
	m_coneArray = QOpenGLBuffer(QOpenGLBuffer::VertexBuffer);
	m_coneArray.create();
	m_coneArray.setUsagePattern(QOpenGLBuffer::StaticDraw);
	m_coneArray.bind();
	m_coneArray.allocate(m_cone.m_vertices.data(), m_cone.m_vertices.size() * sizeof(QVector3D));
	m_coneIndex = QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);
	m_coneIndex.create();
	m_coneIndex.setUsagePattern(QOpenGLBuffer::StaticDraw);
	m_coneIndex.bind();
	m_coneIndex.allocate(m_cone.m_indices.data(), m_cone.m_indices.size() * sizeof(GLushort));
	m_shaderProgram->enableAttributeArray(0);
	m_shaderProgram->setAttributeBuffer(0, GL_FLOAT, 0, 3, sizeof(QVector3D));
	m_coneArray.release();
	m_coneArray.release();
	coneVao->release();
}

void GLcanvas::initArrow()
{
	arrowVao = new QOpenGLVertexArrayObject(this);
	arrowVao->create();
	arrowVao->bind();
	m_arrowArray = QOpenGLBuffer(QOpenGLBuffer::VertexBuffer);
	m_arrowArray.create();
	m_arrowArray.setUsagePattern(QOpenGLBuffer::StaticDraw);
	m_arrowArray.bind();
	m_arrowArray.allocate(m_arrow.m_vertices.data(), m_arrow.m_vertices.size() * sizeof(QVector3D));
	m_arrowIndex = QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);
	m_arrowIndex.create();
	m_arrowIndex.setUsagePattern(QOpenGLBuffer::StaticDraw);
	m_arrowIndex.bind();
	m_arrowIndex.allocate(m_arrow.m_indices.data(), m_arrow.m_indices.size() * sizeof(GLushort));
	
	m_shaderProgram->enableAttributeArray(0);
	m_shaderProgram->setAttributeBuffer(0, GL_FLOAT, 0, 3, 0);
	m_shaderProgram->enableAttributeArray(1);
	m_shaderProgram->setAttributeBuffer(1, GL_FLOAT, m_arrow.m_indices.size() * sizeof(QVector3D), 3, 0);
	
	m_arrowArray.release();
	m_arrowArray.release();
	arrowVao->release();
}

void GLcanvas::panEyeWindow(QVector2D panV)
{
	m_left = m_left + panV.x();
	m_right = m_right + panV.x();
	m_top = m_top + panV.y();
	m_bottom = m_bottom + panV.y();
	m_projectionMatrix.setToIdentity();
	m_projectionMatrix.frustum(m_left, m_right, m_bottom, m_top, m_front, m_back);
	update();
}

void GLcanvas::scaleWordWindow(qreal _scaleFac)
{
	double vpr;          // viewport distortion ratio
	double cx, cy;    // window center
	double sizex, sizey; // window sizes
	double c;

	//Compute canvas viewport ratio.
	vpr = (double)height() / (double) width();


	cx = (m_right + m_left) / 2;
	cy = (m_top + m_bottom) / 2;

	sizex = _scaleFac * abs(m_right - m_left);
	sizey = _scaleFac * abs(m_top - m_bottom);


	c = sizey / sizex;

	if ((c > vpr)) // x é o menor lado logo não é alterado
	{
		sizex = sizey / vpr;
	}
	else if ((c > vpr)) // y é o menor lado logo não é alterado
	{
		sizey = (vpr)*sizey;
	}
	c = (sizey / sizex);

	m_right = cx + (sizex / 2);
	m_left = cx - (sizex / 2);
	m_bottom = cy - (sizey / 2);
	m_top = cy + (sizey / 2);

	m_projectionMatrix.setToIdentity();
	m_projectionMatrix.frustum(m_left, m_right, m_bottom, m_top, m_front, m_back);

	update();
}

void GLcanvas::mousePressEvent(QMouseEvent* event)
{
	m_buttonPressed = true;
	m_pt0 = m_pt1 = QVector2D(event->pos());
	m_mouseButton = event->button();
	QVector2D pt0_eye = convertRasterPtToEye(m_pt0);
	QVector2D pt0_eye_NDC = QVector2D(pt0_eye.x() / m_rightS, pt0_eye.y() / m_topS);
	QMatrix4x4 transform;

	if (m_mouseButton == Qt::LeftButton)
	{
		// build view-projection matrix to compare with screen coordinates
		transform.setToIdentity();
		transform.frustum(m_leftS, m_rightS, m_bottomS, m_topS, m_frontS, m_backS);
		transform.lookAt(m_eye, m_ref, m_up);

		switch (m_currentAction)
		{
		case ActionType::UNDEF:
			break;
		case ActionType::SELECTION:
			m_handleSO->selectObjects(pt0_eye_NDC, &transform);
			if (!m_handleSO->isAnySelected())
			{
				//if no object is selected reset color
				m_handleSO->reset();
				update();
			}

			else
			{
				//if no object is selected reset color
				m_handleSO->changeColorSO();
				update();
			}
			break;
		case ActionType::COLLECTION:
			if (!(m_collector.m_startCollecting))
				m_collector.startCollection();
			switch (m_collector.getObjectType())
			{
			case ObjectType::HALFEDGE:
				if (m_collector.isPointSelected(pt0_eye_NDC, &transform, m_model) )
				{
					//update();
				}
				break;
			case ObjectType::FACE:
				if (m_collector.isEdgeSelected(pt0_eye_NDC, &transform, m_model))
				{
					//update
				}
				break;
			}

			break;
		}
	}

}

void GLcanvas::mouseMoveEvent(QMouseEvent* event)
{
	QVector2D pt1_prev = m_pt1;
	m_pt1 = QVector2D(event->pos());
	// se nenhum botão for pressionado, faça nada
	if (!m_buttonPressed)
	{
		return;
	}
	QVector2D d =  m_pt1 - pt1_prev;
	if (d.length() < m_mouseMoveTol)
	{
		return;
	}
	QVector2D pt1_prev_eye = convertRasterPtToEye(pt1_prev);
	QVector2D pt1_eye = convertRasterPtToEye(m_pt1);
    // treat evemt
	if (m_mouseButton == Qt::MidButton)
	{
		qreal panX = pt1_eye.x() - pt1_prev_eye.x();
		qreal panY = pt1_eye.y() - pt1_prev_eye.y();
		panEyeWindow(QVector2D(-panX, -panY));

	}
	else if (m_mouseButton == Qt::LeftButton)
	{
		qreal mouseMoveFacX = (qreal)(m_pt1.x() - pt1_prev.x()) / (qreal) width();
		qreal mouseMoveFacY = (qreal)(m_pt1.y() - pt1_prev.y()) / (qreal)height();
		if (abs(mouseMoveFacY) > abs(mouseMoveFacX))
		{
			rotateViewAboutX(-mouseMoveFacY);
		}
		else
		{
			rotateViewAboutY(-mouseMoveFacX);

		}

	}
	else if (m_mouseButton == Qt::RightButton)
	{
		QMatrix4x4 projection, camera, rotation;
		QVector2D v1, v2;
		camera.lookAt(m_eye, m_ref, m_up);
		projection.frustum(m_leftS, m_rightS, m_bottomS, m_topS, m_frontS, m_backS);
		QVector2D origin = QVector2D (projection * camera * QVector4D(m_ref, 1));
		v1 = pt1_eye - origin;
		v2 = pt1_prev_eye - origin;
		QVector2D v2Ortho = QVector2D(-v2.y(), v2.x());
		float dot = QVector2D::dotProduct(v1.normalized(), v2.normalized());
		float angle = std::acos(dot) * 180 / M_PI;
		if (QVector2D::dotProduct(v2Ortho, v1) < 0)
		{
			angle = -angle;
		}
		rotation.rotate(-angle, m_zEye);
		m_up = rotation.mapVector(m_up);
		buildEyeSystem();
		update();
	}
	

}

void GLcanvas::mouseReleaseEvent(QMouseEvent* event)
{
	m_buttonPressed = false;
	if (m_currentAction == ActionType::COLLECTION)
	{
		if (m_collector.m_startCollecting)
		{
			switch (m_collector.getObjectType())
			{
			case ObjectType::HALFEDGE:
				if (m_collector.isHedComplete())
				{
					HalfEdge* hed = m_collector.getCollectedHed();
					// reset a cor dos pontos
					hed->getP1()->setColor(Vertex::normalColor);
					hed->getP2()->setColor(Vertex::normalColor);
					m_model->insertNewHed(hed);
					m_collector.reset();
				}
				break;
			case ObjectType::FACE:
				if (m_collector.isFaceComplete())
				{
					Face* face = m_collector.getCollectedFace();
					checkFaceOrientation(face);
					face->reset();
					m_model->insertNewFace(face);
					m_collector.reset();
				}
			}
			update();
		}
	}
}

void GLcanvas::wheelEvent(QWheelEvent* event)
{
	QPoint numDegrees = event->angleDelta();
	if (numDegrees.y() > 0)
	{
		scaleWordWindow(0.9);
	}
	else
	{
		scaleWordWindow(1.1);

	}
	

}

void GLcanvas::keyPressEvent(QKeyEvent* event)
{
	if (event->key() == Qt::Key_H)
	{
		if ((m_currentAction == ActionType::SELECTION) && m_handleSO->m_isFaceSelected)
		{
			for (Face* face : m_handleSO->getSelectedFaces())
			{
				std::vector<Face*>::iterator it = std::find(m_model->m_faces.begin(), m_model->m_faces.end(), face);
				if (it != m_model->m_faces.end())
				{
					m_hiddenFaces.push_back(face);
					m_model->m_faces.erase(it);
				}
				  
			}
			
		}
	}
	else if (event->key() == Qt::Key_R)
	{
		if (m_currentAction == ActionType::SELECTION)
		{
			if (!m_hiddenFaces.empty())
			{
				int nFaces = m_staticFaces.size();
				m_model->m_faces.resize(nFaces);
				for (int i = 0; i < nFaces; i++)
				{
					m_model->m_faces[i] = m_staticFaces[i];
				}
				m_hiddenFaces.clear();
				m_handleSO->m_isFaceSelected = true;

			}
		}
	}
}

void GLcanvas::keyReleaseEvent(QKeyEvent* event)
{
	update();
}

QVector2D GLcanvas::convertRasterPtToEye(QVector2D rasterPt)
{
	
	QVector2D origin(m_left, m_bottom);
	// create scale matrix
	qreal size_w = m_right - m_left;
	qreal size_h = m_top - m_bottom;
	qreal s_x = size_w / qreal(width());
	qreal s_y = size_h / qreal(height());
	QVector2D pt (s_x * rasterPt.x() , s_y * (height() - rasterPt.y() ) );

	return QVector2D(origin + pt);
}

void GLcanvas::rotateViewAboutX(qreal mouseMoveFacY)
{
	qreal angle = mouseMoveFacY * 180.0;
	QMatrix4x4 transform = rotateAboutAfixedPoint(angle, m_xEye, m_ref);
	m_eye = QVector3D (transform.map(QVector4D(m_eye, 1)));
	transform.setToIdentity();
	QQuaternion rot = QQuaternion::fromAxisAndAngle(m_xEye, angle);
	transform.rotate(angle, m_xEye);

	m_up = QVector3D (transform.map(QVector4D(m_up, 0)));
	buildEyeSystem();
	update();
}

void GLcanvas::rotateViewAboutY(qreal mouseMoveFacX)
{
	qreal angle = mouseMoveFacX * 180.0;
	QMatrix4x4 transform = rotateAboutAfixedPoint(angle, m_yEye, m_ref);
	m_eye = QVector3D(transform.map(QVector4D(m_eye, 1)));
	transform.setToIdentity();
	transform.rotate(angle, m_yEye);
	m_up = transform.mapVector(m_up);
	buildEyeSystem();
	update();
}

QMatrix4x4 GLcanvas::rotateAboutAfixedPoint(qreal angle, QVector3D axis, QVector3D point)
{
	QMatrix4x4 transformation;
	transformation.setToIdentity();
	transformation.translate(-point);
	
	//QQuaternion rot = QQuaternion::fromAxisAndAngle(axis, angle);
	transformation.rotate(angle, axis);

	transformation.translate(point);
	return transformation; 
}

void GLcanvas::buildEyeSystem()
{

	m_zEye = -(m_ref - m_eye).normalized();
	m_xEye = QVector3D::crossProduct(m_up, m_zEye).normalized();
	m_yEye = QVector3D::crossProduct(m_zEye, m_xEye);

}

void GLcanvas::checkFaceOrientation(Face* face)
{
	QMatrix4x4 transform;
	QVector3D p1, p2, p3, cross;
	QVector2D p1Clipping, p2Clipping, p3Clipping;
	transform.setToIdentity();
	transform.frustum(m_leftS, m_rightS, m_bottomS, m_topS, m_frontS, m_backS);
	transform.lookAt(m_eye, m_ref, m_up);
	p1 = *face->getHeds()[0]->getP1()->getP();
	p2 = *face->getHeds()[0]->getP2()->getP();
	p3 = * face->getHeds()[1]->getP2()->getP();
	p1Clipping = QVector2D(transform * QVector4D(p1, 1.0));
	p2Clipping = QVector2D(transform * QVector4D(p2, 1.0));
	p3Clipping = QVector2D(transform * QVector4D(p3, 1.0));
	qreal orient = Face::orientedArea(p1Clipping, p2Clipping, p3Clipping);
	if (orient > 0)
	{
		return;
	}
	else 
	{
		face->reverseHeds();
		return;
	}

}

void GLcanvas::fit2world()
{
	qreal viewWindowW, viewWindowH, ratio;
	viewWindowW = viewWindowH = m_max_size;
	ratio = (qreal)height() / (qreal)width();
	if (ratio < 1.0)
	{
		viewWindowW = (1.0 / ratio) * viewWindowH;
	}
	else
	{
		viewWindowH = ratio * viewWindowW;
	}
	m_left = m_leftS =  -0.5 * viewWindowW;
	m_right = m_rightS = 0.5 * viewWindowW;
	m_bottom = m_bottomS =  -0.5 * viewWindowH;
	m_top = m_topS = 0.5 * viewWindowH;
	m_front = m_frontS  = 3.0 * m_max_size;
	m_back = m_backS = 5.0 * m_max_size;
	
	

	// camera
	m_ref = *m_model->getCentroide(); //center
	m_eye = QVector3D(m_ref.x(), m_ref.y(), m_ref.z() + 4.0 * m_max_size);
	m_up = QVector3D(0.0, 1.0, 0.0);

	buildEyeSystem();

	// update
	m_projectionMatrix.setToIdentity();
	m_projectionMatrix.frustum(m_left, m_right, m_bottom, m_top, m_front, m_back);

	update();

}

void GLcanvas::setMouseAction(ActionType actionType)
{
	if (m_currentAction == actionType)
		return;
	m_currentAction = actionType;
}

void GLcanvas::setNewEyeReference(QVector3D* ref)
{
	m_ref = *ref;
	m_eye = QVector3D(m_ref.x(), m_ref.y(), m_ref.z() + 4.0 * m_max_size);
	m_up = QVector3D(0.0, 1.0, 0.0);
}

void GLcanvas::deleteSelectedObjects()
{
	m_handleSO->deleteSelectedObjects();
	update();
}

void GLcanvas::reset()
{
	m_handleSO->reset();
	update();
}




void GLcanvas::setModel(Model* model)
{
	m_model = model;
	

	// create a copy for hidden faces
	
	for (Face* face:m_model->m_faces)
	{
		m_staticFaces.push_back(face);
	}
}

void GLcanvas::setLogTables(QTableWidget* logPointTable, QTableWidget* logEdgeTable, QTableWidget* logFaceTable)
{
	m_logPoints = logPointTable;
	m_logEdges = logEdgeTable;
	m_logFaces = logFaceTable;
	m_handleSO = new HandleSO(m_model, m_logPoints, m_logEdges, m_logFaces);
	QObject::connect(m_logPoints, &QTableWidget::itemChanged, this, &GLcanvas::on_logPoints_itemChanged);
	QObject::connect(m_logPoints, &QTableWidget::itemDoubleClicked, this, &GLcanvas::on_logPoints_doubleClicked);
	QObject::connect(m_logEdges, &QTableWidget::itemChanged, this, &GLcanvas::on_logEdges_itemChanged);
	QObject::connect(m_logEdges, &QTableWidget::itemDoubleClicked, this, &GLcanvas::on_logEdges_doubleClicked);
	QObject::connect(m_logFaces, &QTableWidget::itemChanged, this, &GLcanvas::on_logFaces_itemChanged);
	QObject::connect(m_logFaces, &QTableWidget::itemDoubleClicked, this, &GLcanvas::on_logFaces_doubleClicked);

}
void GLcanvas::on_logPoints_itemChanged(QTableWidgetItem* item)
{
	if (m_handleSO->m_logPointsInit)
	{
		int m = item->row();
		int n = item->column();
		if (!(n == 0))
		{
			QString str = (m_logPoints->item(m, 0)->data(0).toString())[1];
			int id = str.toInt();
			Vertex* v = m_model->m_vertexes[id];
			bool* ok = new bool;
			switch (m_handleSO->getSelectionType())
			{
			case SelectionType::SELECTION:
			{
				qreal value = m_logPoints->item(m, n)->data(0).toReal(ok);
				if (!ok)
				{
					QVector3D p = *v->getP();
					value = p[n - 1];
					m_logPoints->item(m, n)->setData(0, QString::number(value));
				}
				v->setPosValue(n - 1, value);
				break;
			}
			case SelectionType::SELECTIONBC:
			{
				if (n == 1)
				{
					bool bc = v->getBC();
					bool newBC = m_logPoints->item(m, n)->data(0).toBool();
					if (bc != newBC)
					{
						v->setBC(newBC);
						if (newBC)
						{
							v->setTemperature(0);
						}
					}
				}
				else
				{
					if (v->getBC())
					{
						qreal value = m_logPoints->item(m, n)->data(0).toReal(ok);
						if (!ok)
						{
							value = v->getTemperature();
							m_logPoints->item(m, n)->setData(0, QString::number(value));
						}
						v->setTemperature(value);
					}
					break;
				}
			}
			}
			m_handleSO->m_logPointsInit = false;
			delete ok;
			update();
		}
	}
}
void GLcanvas::on_logEdges_itemChanged(QTableWidgetItem* item)
{
	if (m_handleSO->m_logEdgesInit)
	{
		int m = item->row();
		int n = item->column();
		if (!(n == 0))
		{
			QString str = (m_logEdges->item(m, 0)->data(0).toString())[1]; // getId aresta
			int id = str.toInt();
			Edge* e = m_model->m_edges[id];
			bool* ok = new bool;
			unsigned int value = m_logEdges->item(m, n)->data(0).toUInt(ok);
			// se o usario digitar uma string pegue o valor anterior
			if ( (!*ok) || value >= m_model->m_vertexes.size())
			{
				value = e->getIdP1();
				if (n == 2)
				{
					value = e->getIdP2();
				}
				m_logEdges->item(m, n)->setData(0, QString::number(value));
			}
			
			Vertex* v = m_model->m_vertexes[value];
			e->setNewPoint(n, v);
			m_handleSO->m_logEdgesInit = false;
			delete ok;
			update();
		}
	}

}
void GLcanvas::on_logPoints_doubleClicked(QTableWidgetItem* item)
{
	m_handleSO->m_logPointsInit = true;
}

void GLcanvas::on_logEdges_doubleClicked(QTableWidgetItem* item)
{
	m_handleSO->m_logEdgesInit = true;

}

void GLcanvas::on_logFaces_itemChanged(QTableWidgetItem* item)
{
	if (m_handleSO->getSelectionType() == SelectionType::SELECTIONBC)
	{
		if (m_handleSO->m_logFacesInit)
		{
			int m = item->row();
			int n = item->column();
			if (!(n == 0))
			{
				QString str = (m_logFaces->item(m, 0)->data(0).toString())[1];
				int id = str.toInt();
				Face* v = m_model->m_faces[id];
				bool* ok = new bool;
				if (n == 1)
				{
					bool bc = v->getBC();
					bool newBC = m_logFaces->item(m, n)->data(0).toBool();
					if (bc != newBC)
					{
						v->setBC(newBC);
						if (newBC)
						{
							v->setGradient(0);
						}
					}
				}
				else
				{
					if (v->getBC())
					{
						qreal value = m_logFaces->item(m, n)->data(0).toReal(ok);
						if (!ok)
						{
							value = v->getGradient();
							m_logFaces->item(m, n)->setData(0, QString::number(value));
						}
						v->setGradient(value);
					}
				}
				m_handleSO->m_logFacesInit = false;
				delete ok;
				update();
			}

		}
	}

}

void GLcanvas::on_logFaces_doubleClicked(QTableWidgetItem* item)
{
	m_handleSO->m_logFacesInit = true;
}


