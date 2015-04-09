include($$[STARLAB])
include($$[SURFACEMESH])
StarlabTemplate(none)

TEMPLATE = lib
CONFIG += staticlib
QT += opengl

# Build flag
CONFIG(debug, debug|release) {
    CFG = debug
} else {
    CFG = release
}

# Library name and destination
TARGET = HACD
DESTDIR = $$PWD/$$CFG/lib


INCLUDEPATH += "include/"
INCLUDEPATH += "public/"

SOURCES += hacdlib.cpp \
    src/AutoGeometry.cpp \
    src/ConvexDecomposition.cpp \
    src/ConvexHull.cpp \
    src/dgAABBPolygonSoup.cpp \
    src/dgConvexHull3d.cpp \
    src/dgGoogol.cpp \
    src/dgIntersections.cpp \
    src/dgMatrix.cpp \
    src/dgMeshEffect.cpp \
    src/dgPolygonSoupBuilder.cpp \
    src/dgPolyhedra.cpp \
    src/dgQuaternion.cpp \
    src/dgSmallDeterminant.cpp \
    src/dgSphere.cpp \
    src/dgTree.cpp \
    src/dgTypes.cpp \
    src/HACD.cpp \
    src/MergeHulls.cpp \
    src/WuQuantizer.cpp \
    TestHACD.cpp \
    wavefront.cpp

HEADERS += hacdlib.h \
    public/ConvexHull.h \
    public/HACD.h \
    public/JobSwarm.h \
    public/MergeHulls.h \
    public/PlatformConfigHACD.h \
    public/PxVector.h \
    public/SparseArray.h \
    public/WuQuantizer.h \
    wavefront.h

