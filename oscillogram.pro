#-------------------------------------------------
#
# Project created by QtCreator 2023-10-11T15:31:22
#
#-------------------------------------------------

QT       += core gui printsupport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = oscillogram
TEMPLATE = app
DEFINES += QT_DEPRECATED_WARNINGS

CONFIG += c++11

SOURCES += \
        src/main.cpp \
        src/mainwindow.cpp \
        lib/qcustomplot.cpp

HEADERS += \
        src/mainwindow.h \
        lib/qcustomplot.h

FORMS += \
        UI/mainwindow.ui

qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
