TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += $PWD

# FLAGS

SOURCES += main.cpp \
    core/RandomDistributions.cpp \
    core/CODeMDistribution.cpp \
    core/CODeMOperators.cpp \
    core/CODeMProblems.cpp \
    core/UncertaintyKernel.cpp \
    core/utils/LinearInterpolator.cpp \
    libs/DTLZ/DTLZProblems.cpp \
    libs/WFG/ExampleProblems.cpp \
    libs/WFG/ExampleShapes.cpp \
    libs/WFG/ExampleTransitions.cpp \
    libs/WFG/FrameworkFunctions.cpp \
    libs/WFG/Misc.cpp \
    libs/WFG/ShapeFunctions.cpp \
    libs/WFG/TransFunctions.cpp

HEADERS += \
    core/RandomDistributions.h \
    core/CODeMDistribution.h \
    core/CODeMOperators.h \
    core/CODeMProblems.h \
    core/UncertaintyKernel.h \
    core/utils/LinearInterpolator.h \
    libs/DTLZ/DTLZProblems.h \
    libs/WFG/ExampleProblems.h \
    libs/WFG/ExampleShapes.h \
    libs/WFG/ExampleTransitions.h \
    libs/WFG/FrameworkFunctions.h \
    libs/WFG/Misc.h \
    libs/WFG/ShapeFunctions.h \
    libs/WFG/TransFunctions.h
