#include "labeler.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    Labeler w;
    w.show();

    return a.exec();
}
