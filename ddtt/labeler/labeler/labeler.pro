include($$[STARLAB])
include($$[SURFACEMESH])
include($$[OCTREE])
StarlabTemplate(appbundle)

# Build flag for the static libraries
CONFIG(debug, debug|release) {
    CFG = debug
} else {
    CFG = release
}

# NURBS library
LIBS += -L$$PWD/../../NURBS/$$CFG/lib -lNURBS
INCLUDEPATH += ../../NURBS

# Surface Reconstruction library
LIBS += -L$$PWD/../../Reconstruction/$$CFG/lib -lReconstruction
INCLUDEPATH += ../../Reconstruction

# TopoBlender library
LIBS += -L$$PWD/../../StructureGraphLib/$$CFG/lib -lStructureGraphLib
INCLUDEPATH += ../../StructureGraphLib

# Splat Rendering library
LIBS += -L$$PWD/../../GlSplatRendererLib/$$CFG/lib -lGlSplatRendererLib
INCLUDEPATH += ../../GlSplatRendererLib

QT += core gui opengl svg network

TARGET = labeler

FORMS += \
    labeler.ui

HEADERS += \
    labeler.h \
    mydrawarea.h

SOURCES += \
    labeler.cpp \
    main.cpp \
    mydrawarea.cpp
