######################################################################
# Automatically generated by qmake (3.0) Sat Jan 3 11:28:00 2015
######################################################################

TEMPLATE = app
TARGET = scalarfield
INCLUDEPATH += .
QMAKE_CXXFLAGS += -std=c++0x
OBJECTS_DIR = obj

LIBS += -lnoise -L./lib

# Input
SOURCES += src/main.cpp

HEADERS += \
    src/catmulrom.h \
    src/scalarfield.hpp

OTHER_FILES += scalarfied.dox

DISTFILES += $${OTHER_FILES}
