
#include <QtWidgets/QMainWindow>
#include "ui_testeqtmoderngl.h"
#include "model/model.h"
#include "topoSolver/Solid.h"

class testeQtModernGL : public QMainWindow
{
    Q_OBJECT

public:
    testeQtModernGL(QWidget *parent = Q_NULLPTR);

private:
    Ui::testeQtModernGLClass ui;
    Model* m_model;
    topoSolver::Solid* m_solid;
    void onActionFit2World();
    void onActionZoomIn();
    void onActionZoomOut();
    void onActionSelect();
    void onActionPoint();
    void onActionEdge();
    void onActionFace();
    void onActionBoundaryCondition();
    void onActionSolver();
    void onActionDelete();

    void createBoundaryTable();
    void resetTable();
    void createSolid(const int nL);
    void onReadFile();
    void onConvexHull();
};
