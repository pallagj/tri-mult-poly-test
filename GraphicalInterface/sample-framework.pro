# -*- mode: Makefile -*-

TARGET = sample-framework
CONFIG += c++14 qt opengl debug
QT += gui widgets opengl xml
equals (QT_MAJOR_VERSION, 6) {
    QT += openglwidgets
}

HEADERS = MyWindow.h MyViewer.h MyViewer.hpp \
    gorbe-algoritmusok.h \
    bsplinecurve.h \
    vector.h \
    coonssurface.h \
    curve.h \
    External/DMWT/Algorithm/DMWT.h \
    External/DMWT/DataStructure/Boundary.h \
    External/DMWT/DataStructure/CurveInfo.h \
    External/DMWT/DataStructure/EdgeInfo.h \
    External/DMWT/DataStructure/Hole.h \
    External/DMWT/DataStructure/SubKey.h \
    External/DMWT/DataStructure/Tile.h \
    External/DMWT/DataStructure/TriangleInfo.h \
    External/DMWT/Utility/Configure.h \
    External/DMWT/Utility/Point3.h \
    External/DMWT/Utility/Vector3.h \
    External/ExternalLibs/tetgen/tetgen1.4.3/inc/tetgen.h \
    External/ExternalLibs/tetgen/tetgen1.4.3/src/tetgen.h \
    bsplinecurve.h \
    coons_copy.h \
    coonssurface.h \
    curve.h \
    gorbe-algoritmusok.h \
    jet-wrapper.h \
    MyViewer.h \
    MyViewer.hpp \
    MyWindow.h \
    vector.h
SOURCES = MyWindow.cpp MyViewer.cpp main.cpp jet-wrapper.cpp \
    gorbe-algoritmusok.cc \
    bsplinecurve.cpp \
    vector.cpp \
    coonssurface.cpp \
    curve.cpp \
    External/DMWT/Algorithm/DMWT.cpp \
    External/DMWT/DataStructure/Boundary.cpp \
    External/DMWT/DataStructure/CurveInfo.cpp \
    External/DMWT/DataStructure/EdgeInfo.cpp \
    External/DMWT/DataStructure/Hole.cpp \
    External/DMWT/DataStructure/Tile.cpp \
    External/DMWT/DataStructure/TriangleInfo.cpp \
    External/DMWT/Utility/Point3.cpp \
    External/DMWT/Utility/Vector3.cpp \
    External/ExternalLibs/tetgen/tetgen1.4.3/src/predicates.cxx \
    External/ExternalLibs/tetgen/tetgen1.4.3/src/tetgen.cxx \
    External/ExternalLibs/tetgen/tetgen1.4.3/test/example_tetcall.cxx \
    External/ExternalLibs/tetgen/tetgen1.4.3/test/tetcall.cxx \
    bsplinecurve.cpp \
    coonssurface.cpp \
    curve.cpp \
    gorbe-algoritmusok.cc \
    gorbe-kiertekeles.cc \
    jet-wrapper.cpp \
    main.cpp \
    MyViewer.cpp \
    MyWindow.cpp \
    vector.cpp

QMAKE_CXXFLAGS += -O3

unix:INCLUDEPATH += /usr/include/eigen3
unix:LIBS *= -lQGLViewer-qt5 -lOpenMeshCore -lGL -lGLU

RESOURCES = sample-framework.qrc

# Optional
# DEFINES += BETTER_MEAN_CURVATURE
# DEFINES += USE_JET_FITTING
# LIBS += -lCGAL # this library will be header-only from version 5

###########################
# WIN32 instructions only #
###########################

win32 {
    # Replace this variable to the install path of the OpenMesh lib install
    OPENMESH_INSTALL_PATH = 'C:\Program Files\OpenMesh'

    # Replace this variable to the install path of libQGLViewer
    LIBQGLVIEWER_INSTALL_PATH = 'C:\Program Files\libQGLViewer'

    # If your OpenMesh source is separate from the lib install replace this variable
    OPENMESH_SRC_INSTALL_PATH = $$OPENMESH_INSTALL_PATH\include\

    DEFINES += NOMINMAX _USE_MATH_DEFINES

    INCLUDEPATH += '$$OPENMESH_SRC_INSTALL_PATH' '$$LIBQGLVIEWER_INSTALL_PATH'

    LIBS += -lOpenGL32 -lGLU32
    LIBS += -L'$$OPENMESH_INSTALL_PATH\lib' -L'$$LIBQGLVIEWER_INSTALL_PATH\QGLViewer'
    Release:LIBS += -lOpenMeshCore -lQGLViewer2
    Debug:LIBS += -lOpenMeshCored -lQGLViewerd2
}
