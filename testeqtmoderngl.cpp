#include <QMessageBox>
#include <QThread>
#include <QProgressDialog>
#include <QFileDialog>
#include "testeqtmoderngl.h"
#include "dialogs/insertpointdialog.h"
#include "dialogs/solverDialog.h"
#include "topoSolver/Edge.h"
#include "topoSolver/Vertex.h"
#include "topoSolver/Face.h"
#include "topoSolver/HalfEdge.h"
#include "fmm/FMM.h"
#include "geom/convexHull.h"

testeQtModernGL::testeQtModernGL(QWidget *parent)
    : QMainWindow(parent)
{
    ui.setupUi(this);
    m_model = new Model();
    ui.glDisplay->setModel(m_model);
    ui.glDisplay->setLogTables(ui.logPoints, ui.logEdges, ui.logFaces);
    //
		ui.actionSelect->setChecked(true);
		ui.actionPoint->setChecked(false);
		ui.actionEdge->setChecked(false);
		ui.actionFace->setChecked(false);
		//
	connect(ui.actionFit2World, &QAction::triggered, this, &testeQtModernGL::onActionFit2World);
	connect(ui.actionZoom_in, &QAction::triggered, this, &testeQtModernGL::onActionZoomIn);
	connect(ui.actionZoom_out, &QAction::triggered, this, &testeQtModernGL::onActionZoomOut);
	connect(ui.actionSelect, &QAction::triggered, this, &testeQtModernGL::onActionSelect);
	connect(ui.actionPoint, &QAction::triggered, this, &testeQtModernGL::onActionPoint);
	connect(ui.actionEdge, &QAction::triggered, this, &testeQtModernGL::onActionEdge);
	connect(ui.actionFace, &QAction::triggered, this, &testeQtModernGL::onActionFace);
	connect(ui.actionBoundaryCondition, &QAction::triggered, this, &testeQtModernGL::onActionBoundaryCondition);
	connect(ui.actionSolver, &QAction::triggered, this, &testeQtModernGL::onActionSolver);
	connect(ui.actionDelete, &QAction::triggered, this, &testeQtModernGL::onActionDelete);
	connect(ui.actionOpen, &QAction::triggered, this, &testeQtModernGL::onReadFile);
	connect(ui.actionConvexHull, &QAction::triggered, this, &testeQtModernGL::onConvexHull);



}

void testeQtModernGL::onActionFit2World()
{
	ui.glDisplay->fit2world();
}

void testeQtModernGL::onActionZoomIn()
{
	ui.glDisplay->scaleWordWindow(0.9);
}

void testeQtModernGL::onActionZoomOut()
{
	ui.glDisplay->scaleWordWindow(1.1);
}

void testeQtModernGL::onActionSelect()
{
	ui.actionSelect->setChecked(true);
	ui.actionEdge->setEnabled(true);
	ui.actionFace->setEnabled(true);
	ui.actionEdge->setChecked(false);
	ui.actionFace->setChecked(false);
	ui.actionDelete->setEnabled(true);
	ui.glDisplay->setMouseAction(ActionType::SELECTION);
	ui.glDisplay->m_handleSO->setSelectionType(SelectionType::SELECTION);
	if (ui.actionBoundaryCondition->isChecked())
	{
		resetTable();
		ui.actionBoundaryCondition->setChecked(false);
	}
}

void testeQtModernGL::onActionPoint()
{
	ui.actionSelect->setChecked(true);
	ui.actionEdge->setChecked(false);
	ui.actionFace->setChecked(false);
	ui.glDisplay->setMouseAction(ActionType::SELECTION);

	InsertPointDialog dg{ this };
	if (dg.exec() == QDialog::Accepted)
	{
		QVector3D pt = dg.getPoint();
		m_model->insertNewVertex(&pt);
		ui.glDisplay->setNewEyeReference(m_model->getCentroide());
	}
	ui.glDisplay->update();
}

void testeQtModernGL::onActionEdge()
{
	ui.actionSelect->setChecked(false);
	ui.actionEdge->setChecked(true);
	ui.actionFace->setChecked(false);
	ui.actionDelete->setEnabled(false);

	ui.glDisplay->setMouseAction(ActionType::COLLECTION);
	ui.glDisplay->m_collector.setObjectType(ObjectType::HALFEDGE);
	ui.glDisplay->reset();
	ui.glDisplay->m_collector.reset();




}

void testeQtModernGL::onActionFace()
{
	ui.actionSelect->setChecked(false);
	ui.actionEdge->setChecked(false);
	ui.actionFace->setChecked(true);
	ui.actionDelete->setEnabled(false);

	ui.glDisplay->setMouseAction(ActionType::COLLECTION);
	ui.glDisplay->m_collector.setObjectType(ObjectType::FACE);
	ui.glDisplay->reset();
	ui.glDisplay->m_collector.reset();

}

void testeQtModernGL::onActionBoundaryCondition()
{
	int V = m_model->m_vertexes.size();
	int A = m_model->m_edges.size();
	int F = m_model->m_faces.size();
	if (V + F - A == 2 && F != 0)
	{
		ui.actionSelect->setChecked(false);
		ui.actionBoundaryCondition->setChecked(true);
		ui.actionPoint->setEnabled(false);
		ui.actionEdge->setEnabled(false);
		ui.actionFace->setEnabled(false);
		ui.actionDelete->setEnabled(false);
		ui.actionBoundaryCondition->setEnabled(true);
		ui.glDisplay->setMouseAction(ActionType::SELECTION);
		ui.glDisplay->m_handleSO->setSelectionType(SelectionType::SELECTIONBC);
		ui.actionSolver->setEnabled(true);
		createBoundaryTable();
	}	else 
	{
		ui.actionBoundaryCondition->setChecked(false);
		QMessageBox::about(this, "Numero de Euler invalido", " Malha criada invalida. Por favor, verificar numero de Euler");
	}
}



void testeQtModernGL::onActionSolver()
{
	ui.actionBoundaryCondition->setChecked(false);
	SolverDialog dg{ m_model , this  };

	if (dg.exec() == QDialog::Accepted)
	{		
		QProgressDialog progress("calculando mapa de calor...", "Abort ", 0, 100, this);
		progress.setWindowModality(Qt::WindowModal);
		progress.setWindowTitle("Solver");
		progress.show();	
		for (int i = 25; i < 50; i++)
		{
			progress.setValue(i);
			QThread::msleep(50);
		}

		createSolid(dg.getLevel());
		const int NG = 20; // gauss points to integrate fmm
		FMM fmm(m_solid, dg.getTruncamentoTerm(), NG);
		for (int i = 25; i < 70; i++)
		{
			progress.setValue(i);
			QThread::msleep(50);
		}
		//fmm.showX();
	}
}

void testeQtModernGL::onActionDelete()
{
	ui.glDisplay->deleteSelectedObjects();
}

void testeQtModernGL::createBoundaryTable()
{
	QTableWidgetItem* ___qtablewidgetitem1;
	QTableWidgetItem* ___qtablewidgetitem2;

	ui.logPoints->setColumnCount(3);
	___qtablewidgetitem1 = ui.logPoints->horizontalHeaderItem(1);
	___qtablewidgetitem1->setText(QCoreApplication::translate("testeQtModernGLClass", "B.C.", nullptr));
	___qtablewidgetitem2 = ui.logPoints->horizontalHeaderItem(2);
	___qtablewidgetitem2->setText(QCoreApplication::translate("testeQtModernGLClass", "temperatura", nullptr));
	
	ui.logPoints->setMinimumSize(QSize(150, 200));
	ui.logPoints->setMaximumSize(QSize(250, 300));

	ui.logFaces->setColumnCount(3);
	___qtablewidgetitem1 = ui.logFaces->horizontalHeaderItem(1);
	___qtablewidgetitem1->setText(QCoreApplication::translate("testeQtModernGLClass", "B.C.", nullptr));
	___qtablewidgetitem2 = ui.logFaces->horizontalHeaderItem(2);
	___qtablewidgetitem2->setText(QCoreApplication::translate("testeQtModernGLClass", "gradiente", nullptr));
	
	ui.logFaces->setMinimumSize(QSize(125, 200));
	ui.logFaces->setMaximumSize(QSize(200, 300));



	ui.logEdges->setMinimumSize(QSize(150, 0));
	ui.logEdges->setMaximumSize(QSize(250, 0));
}

void testeQtModernGL::resetTable()
{
	ui.logPoints->setColumnCount(4);
	QTableWidgetItem* __qtablewidgetitem = new QTableWidgetItem();
	ui.logPoints->setHorizontalHeaderItem(3, __qtablewidgetitem);
	QTableWidgetItem* ___qtablewidgetitem1 = ui.logPoints->horizontalHeaderItem(1);
	___qtablewidgetitem1->setText(QCoreApplication::translate("testeQtModernGLClass", "x", nullptr));
	QTableWidgetItem* ___qtablewidgetitem2 = ui.logPoints->horizontalHeaderItem(2);
	___qtablewidgetitem2->setText(QCoreApplication::translate("testeQtModernGLClass", "y", nullptr));
	QTableWidgetItem* ___qtablewidgetitem3 = ui.logPoints->horizontalHeaderItem(3);
	___qtablewidgetitem3->setText(QCoreApplication::translate("testeQtModernGLClass", "z", nullptr));

	ui.logFaces->setColumnCount(4);
	QTableWidgetItem* __qtablewidgetitemF = new QTableWidgetItem();
	ui.logFaces->setHorizontalHeaderItem(3, __qtablewidgetitemF);
	QTableWidgetItem* ___qtablewidgetitem8 = ui.logFaces->horizontalHeaderItem(1);
	___qtablewidgetitem8->setText(QCoreApplication::translate("testeQtModernGLClass", "A0", nullptr));
	QTableWidgetItem* ___qtablewidgetitem9 = ui.logFaces->horizontalHeaderItem(2);
	___qtablewidgetitem9->setText(QCoreApplication::translate("testeQtModernGLClass", "A1", nullptr));
	QTableWidgetItem* ___qtablewidgetitem10 = ui.logFaces->horizontalHeaderItem(3);
	___qtablewidgetitem10->setText(QCoreApplication::translate("testeQtModernGLClass", "A2", nullptr));

	ui.logPoints->setMinimumSize(QSize(150, 105));
	ui.logPoints->setMaximumSize(QSize(250, 200));

	ui.logFaces->setMinimumSize(QSize(150, 105));
	ui.logFaces->setMaximumSize(QSize(250, 200));

	ui.logEdges->setMinimumSize(QSize(150, 105));
	ui.logEdges->setMaximumSize(QSize(250, 200));
}

void testeQtModernGL::createSolid(const int nL)
{
	const int nHeds = m_model->m_heds.size();
	const int nVerts = m_model->m_vertexes.size();
	const int nEdges = m_model->m_edges.size();
	const int nEl = m_model->m_faces.size();
	// Initialize
	std::vector<topoSolver::Vertex*> vertexs(nVerts);
	std::vector<topoSolver::Vertex*> vertexsBd;
	std::vector<topoSolver::Edge*> edges(nEdges);
	std::vector<topoSolver::HalfEdge*> heds(nHeds);
	std::vector<topoSolver::Face*> elements(nEl);

	for (topoApp::Edge* edge : m_model->m_edges)
	{
		const int id = edge->getId();
		edges[id] = new topoSolver::Edge(id);
	}

	for (topoApp::Vertex* pt : m_model->m_vertexes)
	{
		const int id = pt->getId();
		QVector3D* coord = pt->getP();
		vertexs[id] = new topoSolver::Vertex(Point(coord->x(), coord->y(), coord->z()), id);
		if (pt->getBC())
		{
			vertexs[id]->setBC((double)pt->getTemperature());
			vertexsBd.push_back(vertexs[id]);
		}
	}

	for (topoApp::HalfEdge* hed : m_model->m_heds)
	{
		int inc[2] = { hed->getP1()->getId(), hed->getP2()->getId() };
		const int id = hed->getId();
		const int edgeId = hed->getEdge()->getId();
		heds[id] = new topoSolver::HalfEdge(inc, id, edges[edgeId], vertexs[inc[0]], vertexs[inc[1]]);
		edges[edgeId]->insertHed(heds[id]);

	}

	for (topoApp::Face* face : m_model->m_faces)
	{
		std::vector<topoApp::HalfEdge*> faceHedsApp = face->get3Heds();
		topoSolver::HalfEdge* faceHeds[3] = { heds[faceHedsApp[0]->getId()] , heds[faceHedsApp[1]->getId()], heds[faceHedsApp[2]->getId()] };
		elements[face->getId()] = new topoSolver::Face( faceHedsApp[0]->getId(), (double) face->getGradient(), faceHeds);
		for (int i = 0; i < 3; i++)
		{
			faceHeds[i]->m_heNext = faceHeds[(i + 1) % 3];
			faceHeds[i]->m_el = elements[face->getId()];
		}
	}
	m_solid = new topoSolver::Solid(vertexs, edges, heds, elements, nL, vertexsBd);



}

void testeQtModernGL::onReadFile()
{
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Open File"), "C:/Users/hilto/OneDrive/Documentos/Estudos de Contorno/Entradas", tr("Image Files (*.txt)"));

	QFile file(fileName);
	if (!file.open(QIODevice::ReadOnly)) {
		QMessageBox::information(0, "error", file.errorString());
	}

	QTextStream in(&file);

	while (!in.atEnd()) {
		QString line = in.readLine();
		QStringList fields = line.split(",");
		QVector3D pt(fields[0].toDouble(), fields[1].toDouble(), fields[2].toDouble());
		m_model->insertNewVertex(&pt);
	}
	file.close();
	ui.glDisplay->setNewEyeReference(m_model->getCentroide());
	ui.glDisplay->update();

	

}

void testeQtModernGL::onConvexHull()
{
	ConvexHull hull(m_model);
	hull.initModel();
	ui.glDisplay->update();
}

