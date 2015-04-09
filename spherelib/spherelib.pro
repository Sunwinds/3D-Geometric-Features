include($$[STARLAB])
include($$[SURFACEMESH])
StarlabTemplate(none)

QT += opengl

TARGET = spherelib
TEMPLATE = lib
CONFIG += staticlib

SOURCES += spherelib.cpp
HEADERS += spherelib.h

# Build options
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}
DESTDIR = $$PWD/$$CFG/lib

# Enable debug in release mode
win32{
    QMAKE_CXXFLAGS_RELEASE += /Zi
    QMAKE_LFLAGS_RELEASE += /DEBUG
}
