include($$[STARLAB])
include($$[SURFACEMESH])
include($$[OCTREE])
include($$[NANOFLANN])
StarlabTemplate(plugin)

QT += gui opengl xml svg

HEADERS += \
    particles.h \
    particles-widget.h \
    ParticleMesh.h \
    Particle.h \
    Raytracing.h

SOURCES += \
    particles.cpp \
    particles-widget.cpp \
    ParticleMesh.cpp \
    Particle.cpp \
    Raytracing.cpp

FORMS       += particles-widget.ui
RESOURCES   += particles.qrc

# Build options
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}

# Sphere library
LIBS += -L$$PWD/../spherelib/$$CFG/lib -lspherelib
INCLUDEPATH += ../spherelib

# External library
win32:LIBS += -L"$$_PRO_FILE_PWD_/embree2/"
