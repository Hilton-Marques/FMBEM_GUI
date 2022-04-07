#include "solverDialog.h"
#include <cmath>
SolverDialog::SolverDialog(Model* model_ ,  QWidget* parent) :
  QDialog(parent), ui(new Ui::solverDialog)
{
  ui->setupUi(this);
  //init values
  m_nL = 4;
  m_N = 3;
  ui->lineEdit->setText("4");
  ui->lineEdit_1->setText("3");
  m_model = model_;
  m_nEl = m_model->m_faces.size();
  m_nVerts = m_model->m_vertexes.size();
  m_nEdges = m_model->m_edges.size();
  QObject::connect(ui->lineEdit, &QLineEdit::editingFinished, this, &SolverDialog::on_lineEdit_changed);
  QObject::connect(ui->lineEdit_1, &QLineEdit::editingFinished, this, &SolverDialog::on_lineEdit_1_changed);

  update();
}

void SolverDialog::setModel(Model* model_)
{
  m_model = model_;
}

void SolverDialog::update()
{
  const int nTotalEl = m_nEl * (std::pow(4, m_nL) - 1) / 3; // sum of P.G.
  const int nTotalVerts = calculateNMaxVerts(m_nL);
  ui->textEdit->setText(QString::number(nTotalVerts));
  ui->textEdit_1->setText(QString::number(nTotalEl));
}

int SolverDialog::calculateNMaxVerts(const int nL)
{
  int m_nL = 2;
  int maxNVerts = m_nVerts + m_nEdges;
  int nEdges = m_nEdges;
  while (m_nL < nL)
  {
    nEdges = nEdges * 2 + m_nEl * std::pow(4, m_nL - 2) * 3;
    maxNVerts = maxNVerts + nEdges;
    m_nL++;
  }
  return maxNVerts;
}
void SolverDialog::on_lineEdit_changed()
{
  bool* ok = new bool;
  double value = ui->lineEdit->text().toInt(ok);
  m_nL = value;
  if (!(*ok))
  {
    ui->lineEdit->setText("4");
    m_nL = 4;
  }
  update();
  delete ok;
}


void SolverDialog::on_lineEdit_1_changed()
{
  bool* ok = new bool;
  double value = ui->lineEdit_1->text().toInt(ok);
  m_N = value;
  if (!(*ok))
  {
    ui->lineEdit_1->setText("3");
    m_N = 3;
  }
  delete ok;
}
