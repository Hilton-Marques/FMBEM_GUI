#include "insertpointdialog.h"


InsertPointDialog::InsertPointDialog(QWidget* parent) :
  QDialog(parent), ui(new Ui::InsertPointDialog)
{
  ui->setupUi(this);
  QObject::connect(ui->lineEdit, &QLineEdit::editingFinished, this, &InsertPointDialog::on_lineEdit1_changed);
}

void InsertPointDialog::on_lineEdit1_changed()
{
  
  bool* ok = new bool; 
  double value = ui->lineEdit->text().toFloat(ok);
  if (!(*ok))
  {
    ui->lineEdit->setText("0.00");
  }
  delete ok;

}

void InsertPointDialog::on_lineEdit2_changed()
{
  bool* ok = new bool;
  double value = ui->lineEdit_2->text().toFloat(ok);
  if (!(*ok))
  {
    ui->lineEdit_2->setText("0.00");
  }
  delete ok;
}

void InsertPointDialog::on_lineEdit3_changed()
{
  bool* ok = new bool;
  double value = ui->lineEdit_3->text().toFloat(ok);
  if (!(*ok))
  {
    ui->lineEdit_3->setText("0.00");
  }
  delete ok;
}


QVector3D InsertPointDialog::getPoint()
{
  float x = ui->lineEdit->text().toFloat();
  float y = ui->lineEdit_2->text().toFloat();
  float z = ui->lineEdit_3->text().toFloat();

  return QVector3D(x,y,z);
}

