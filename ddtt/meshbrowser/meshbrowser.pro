include($$[STARLAB])
include($$[SURFACEMESH])
include($$[OCTREE])
StarlabTemplate(appbundle)

QT += core gui opengl svg network

TARGET = meshbrowser

HEADERS += meshbrowser.h mydrawarea.h
SOURCES += meshbrowser.cpp main.cpp  mydrawarea.cpp
FORMS += meshbrowser.ui

RC_FILE = meshbrowser.rc
