TEMPLATE = lib
CONFIG += staticlib

# Build flag
CONFIG(debug, debug|release) {
    CFG = debug
} else {
    CFG = release
}

# Library name and destination
TARGET = AuctionLib
DESTDIR = $$PWD/$$CFG/lib

SOURCES += \
    Assignment.cpp \
    iAuction.cpp \
    main.cpp \
    Utils.cpp \
    BipartiteGraph.cpp \
    Hungarian.cpp \
    PlotGraph.cpp

HEADERS += \
    Assignment.h \
    CmdParser.h \
    Define.h \
    iAuction.h \
    Utils.h \
    getopt.h \
    BipartiteGraph.h \
    Hungarian.h \
    Matrix.h \
    PlotGraph.h
