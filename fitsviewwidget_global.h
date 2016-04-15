#ifndef FITSVIEWWIDGET_GLOBAL_H
#define FITSVIEWWIDGET_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(FITSVIEWWIDGET_LIBRARY)
#  define FITSVIEWWIDGETSHARED_EXPORT Q_DECL_EXPORT
#else
#  define FITSVIEWWIDGETSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // FITSVIEWWIDGET_GLOBAL_H
