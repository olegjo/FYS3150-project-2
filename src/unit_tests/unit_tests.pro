TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    ../jacobi.cpp \
    ../potentials.cpp \
    unit_tests.cpp

HEADERS += \
    ../jacobi.h \
    ../potentials.h \
    unit_tests.h
