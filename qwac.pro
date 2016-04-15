#-------------------------------------------------
#
# Project created by QtCreator 2016-04-06T10:48:06
#
#-------------------------------------------------

QT       += widgets

TARGET = qwac
TEMPLATE = lib

DEFINES += QWAC_LIBRARY

QMAKE_CXXFLAGS += -std=c++11

SOURCES += \
    FitsViewWidget.cpp \
    focuswidget.cpp \
    expparamsdialog.cpp \
    focussingplotdialog.cpp

HEADERS += \
    FitsViewWidget.h \
    fitsviewwidget_global.h \
    focuswidget_global.h \
    focuswidget.h \
    expparamsdialog.h \
    focussingplotdialog.h


FORMS += \
    focuswidget.ui \
    expParamsForm.ui \
    focussingCurvePlot.ui

unix {
    target.path = /usr/lib
    INSTALLS += target
}


unix:!macx: LIBS += -L/usr/lib64/ -lcfitsio

INCLUDEPATH += /usr/include
DEPENDPATH += /usr/include



unix: CONFIG += link_pkgconfig
unix: PKGCONFIG += Qt5Qwt
