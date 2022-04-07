#include "testeqtmoderngl.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    testeQtModernGL w;
    w.show();
    return a.exec();
}
