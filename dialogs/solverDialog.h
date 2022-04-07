#ifndef SOLVERDIALOG_H
#define SOLVERDIALOG_H

#include "ui_dialogSolver.h"
#include "../model/model.h"
#include <QDialog>

class SolverDialog : public QDialog
{
  Q_OBJECT
public:
  SolverDialog(Model* model , QWidget* parent = Q_NULLPTR);
  void setModel(Model* model_);
  void update();
  int calculateNMaxVerts(const int nL);
  int getLevel() { return m_nL; }; 
  int getTruncamentoTerm() { return m_N;  };
private slots:
  void on_lineEdit_changed();
  void on_lineEdit_1_changed();


private:
  Ui::solverDialog* ui;
  Model* m_model;
  int m_nL; // numero de niveis
  int m_nEl; // numero de faces
  int m_nVerts; // numero de vertices
  int m_nEdges;
  int m_N;
};



#endif // !SOLVERDIALOG_H
