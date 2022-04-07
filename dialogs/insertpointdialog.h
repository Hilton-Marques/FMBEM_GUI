#ifndef INSERT_POINT_DIALOG_H
#define INSERT_POINT_DIALOG_H

#include "ui_DialogButtonBottom.h"
#include <QDialog>
#include <QVector3D>

class InsertPointDialog : public QDialog
{
	Q_OBJECT
public:
	InsertPointDialog(QWidget* parent = Q_NULLPTR);
  
	QVector3D getPoint();

private slots:
	void on_lineEdit1_changed();
	void on_lineEdit2_changed();
	void on_lineEdit3_changed();
private:
	Ui::InsertPointDialog* ui;
};


#endif // !INSERT_POINT_DIALOG_H
