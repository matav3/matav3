#include <QtGui/QApplication>
#include "mainwindow.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    QStringList args = a.arguments();
    MainWindow w;
    w.show();
    if (args.count() == 2){
        std::cerr << "argument required" << endl;
        w.parameterLine(args[1]);
    }
    return a.exec();
}
