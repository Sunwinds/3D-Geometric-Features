include($$[STARLAB])
include($$[SURFACEMESH])
include($$[OCTREE])
include($$[CHOLMOD])

StarlabTemplate(plugin)

# Build flag
CONFIG(debug, debug|release) {
    CFG = debug
} else {
    CFG = release
}

HEADERS += BatchProc.h \
    SegMeshLoader.h \
	GeoHeatHelper.h \
    principal_curvature.h \
    Monge_via_jet_fitting.h \
	Decimater.h \
	SimpleMatrix.h \
    emd.h
SOURCES += BatchProc.cpp \
    SegMeshLoader.cpp \
    principal_curvature.cpp \
    Monge_via_jet_fitting.cpp \
    emd.cpp

#RESOURCES += \
    #BatchProc.qrc

#FORMS += \
    #BatchProcWidget.ui
