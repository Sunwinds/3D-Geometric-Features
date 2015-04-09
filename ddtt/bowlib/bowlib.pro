# Configuration
TEMPLATE = lib
CONFIG += staticlib

# Build flag
CONFIG(debug, debug|release) {
    CFG = debug
} else {
    CFG = release
}

TARGET = bow
DESTDIR = $$PWD/$$CFG/lib

SOURCES +=  galif.cpp \
            vocabulary.cpp \
            BofSearchManager.cpp \
            inverted_index.cpp \
            tf_idf.cpp

HEADERS +=  galif.h \
            types.h \
            image_sampler.h \
            utilities.h \
            vocabulary.h \
            kmeans.h \
            BofSearchManager.h \
            inverted_index.h \
            tf_idf.h \
            quantizer.h \
			histvw.h \
            iofiles.h
			
# OpenCV
win32{
    INCLUDEPATH *= C:/Development/opencv/build/include
    LIBS *= -L"C:/Development/opencv/build/x64/vc11/lib"
    OpenCV_VERSION = 249
    Debug:LIBS *= -lopencv_core$${OpenCV_VERSION}d -lopencv_highgui$${OpenCV_VERSION}d -lopencv_features2d$${OpenCV_VERSION}d -lopencv_nonfree$${OpenCV_VERSION}d -lopencv_flann$${OpenCV_VERSION}d -lopencv_imgproc$${OpenCV_VERSION}d
    Release:LIBS *= -lopencv_core$${OpenCV_VERSION} -lopencv_highgui$${OpenCV_VERSION} -lopencv_features2d$${OpenCV_VERSION} -lopencv_nonfree$${OpenCV_VERSION} -lopencv_flann$${OpenCV_VERSION} -lopencv_imgproc$${OpenCV_VERSION}
}

# OpenMP
win32{
    QMAKE_CXXFLAGS *= /openmp
    QMAKE_CXXFLAGS *= /MP
}
unix:!mac{
    QMAKE_CXXFLAGS *= -fopenmp
    LIBS += -lgomp
}

# Other
win32{
    DEFINES *= NOMINMAX _USE_MATH_DEFINES
}

QMAKE_CXXFLAGS_RELEASE += /Zi
QMAKE_LFLAGS_RELEASE += /DEBUG
