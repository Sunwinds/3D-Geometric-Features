include($$[STARLAB])
include($$[SURFACEMESH])
include($$[CHOLMOD])

StarlabTemplate(plugin)

# Build flag
CONFIG(debug, debug|release) {
    CFG = debug
} else {
    CFG = release
}

# HACD library
LIBS += -L$$PWD/../HACDLib/$$CFG/lib -lHACD
INCLUDEPATH += ../HACDLib

HEADERS += hacd-plugin.h
SOURCES += hacd-plugin.cpp
