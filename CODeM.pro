TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += $PWD

# FLAGS

SOURCES += main.cpp \
    core/distributions/IDistribution.cpp \
    core/distributions/LinearDistribution.cpp \
    core/distributions/MergedDistribution.cpp \
    core/distributions/PeakDistribution.cpp \
    core/distributions/SampledDistribution.cpp \
    core/distributions/UniformDistribution.cpp \
    core/utils/AbstractInterpolator.cpp \
    core/utils/LinearInterpolator.cpp \
    core/CODeM1.cpp \
    core/CODeM2.cpp \
    core/CODeM3.cpp \
    core/CODeM4.cpp \
    core/CODeM5.cpp \
    core/CODeM6.cpp \
    core/CODeMDistribution.cpp \
    core/CODeMOperators.cpp \
    core/CODeMProblems.cpp \
    core/UncertaintyKernel.cpp \
    libs/DTLZ/DTLZProblems.cpp \
    libs/WFG/ExampleProblems.cpp \
    libs/WFG/ExampleShapes.cpp \
    libs/WFG/ExampleTransitions.cpp \
    libs/WFG/FrameworkFunctions.cpp \
    libs/WFG/Misc.cpp \
    libs/WFG/ShapeFunctions.cpp \
    libs/WFG/TransFunctions.cpp

#include(deployment.pri)
#qtcAddDeployment()

HEADERS += \
    core/distributions/IDistribution.h \
    core/distributions/LinearDistribution.h \
    core/distributions/MergedDistribution.h \
    core/distributions/PeakDistribution.h \
    core/distributions/SampledDistribution.h \
    core/distributions/UniformDistribution.h \
    core/utils/AbstractInterpolator.h \
    core/utils/LinearInterpolator.h \
    core/CODeM1.h \
    core/CODeM2.h \
    core/CODeM3.h \
    core/CODeM4.h \
    core/CODeM5.h \
    core/CODeM6.h \
    core/CODeMDistribution.h \
    core/CODeMOperators.h \
    core/CODeMProblems.h \
    core/UncertaintyKernel.h \
    libs/DTLZ/DTLZProblems.h \
    libs/WFG/ExampleProblems.h \
    libs/WFG/ExampleShapes.h \
    libs/WFG/ExampleTransitions.h \
    libs/WFG/FrameworkFunctions.h \
    libs/WFG/Misc.h \
    libs/WFG/ShapeFunctions.h \
    libs/WFG/TransFunctions.h

