#include "imagebrowser.h"
#include <QApplication>

#include <windows.h>

int main(int argc, char *argv[])
{
	SetErrorMode(SEM_NOGPFAULTERRORBOX); 
	SetErrorMode(SEM_FAILCRITICALERRORS | SetErrorMode (0)); 

    QApplication a(argc, argv);

    ImageBrowser w;
    w.showMaximized();

    return a.exec();
}
