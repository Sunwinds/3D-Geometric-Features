include($$[STARLAB])
include($$[SURFACEMESH])
include($$[CHOLMOD])

include($$[NANOFLANN])
include($$[OCTREE])
StarlabTemplate(none)

TEMPLATE = lib
CONFIG += staticlib
QT += opengl xml svg concurrent

# Build flag
CONFIG(debug, debug|release) {
    CFG = debug
} else {
    CFG = release

    DEFINES += QT_NO_DEBUG
    DEFINES += QT_NO_DEBUG_OUTPUT
}

# Library name and destination
TARGET = StructureGraphLib
DESTDIR = $$PWD/$$CFG/lib

# NURBS library
LIBS += -L$$PWD/../NURBS/$$CFG/lib -lNURBS
INCLUDEPATH += ../NURBS

# Surface Reconstruction library
LIBS += -L$$PWD/../Reconstruction/$$CFG/lib -lReconstruction
INCLUDEPATH += ../Reconstruction

# Splat Rendering library
LIBS += -L$$PWD/../GlSplatRendererLib/$$CFG/lib -lGlSplatRendererLib
INCLUDEPATH += ../GlSplatRendererLib

HEADERS += StructureNode.h \
    StructureGraph.h \
    StructureCurve.h \
    StructureSheet.h \
    StructureLink.h \
    GraphEmbed.h \
    TopoBlender.h \
    DynamicGraph.h \
    DynamicGraphGlobal.h \
    GraphDistance.h \
    ExportDynamicGraph.h \
    GraphCorresponder.h \
    Scheduler.h \
    Task.h \
    SchedulerWidget.h \
    TimelineSlider.h \
    StructureGlobal.h \
    Synthesizer.h \
    SynthesisManager.h \
    Sampler.h \
    SimilarSampling.h \
    SpherePackSampling.h \
    AbsoluteOrientation.h \
    TaskCurve.h \
    TaskSheet.h \
    Relink.h \
    GraphModifyWidget.h \
    GraphDissimilarity.h \
    GraphExplorer.h

SOURCES += StructureGraph.cpp \
    StructureCurve.cpp \
    StructureSheet.cpp \
    StructureLink.cpp \
    TopoBlender.cpp \
    DynamicGraph.cpp \
    GraphDistance.cpp \
    GraphCorresponder.cpp \
    Scheduler.cpp \
    Task.cpp \
    SchedulerWidget.cpp \
    TimelineSlider.cpp \
    Synthesizer.cpp \
    SynthesisManager.cpp \
    Sampler.cpp \
    SimilarSampling.cpp \
    AbsoluteOrientation.cpp \
    TaskCurve.cpp \
    TaskSheet.cpp \
    Relink.cpp \
    GraphModifyWidget.cpp \
    GraphDissimilarity.cpp \
    GraphExplorer.cpp

# Graph visualization
SOURCES += QGraphViz/svgview.cpp
HEADERS += QGraphViz/svgview.h

FORMS += SchedulerWidget.ui GraphModifyWidget.ui GraphExplorer.ui

mac:QMAKE_CXXFLAGS += -fopenmp
mac:QMAKE_LFLAGS += -fopenmp

unix:!mac:QMAKE_CXXFLAGS = $$QMAKE_CFLAGS -fpermissive
unix:!mac:LIBS += -lGLU

win32:QMAKE_CXXFLAGS *= /wd4267 /wd4005
win32:QMAKE_CXXFLAGS_RELEASE *= /Zi
#win32:QMAKE_CXXFLAGS_RELEASE *= /Od
win32:QMAKE_LFLAGS_RELEASE *= /DEBUG
