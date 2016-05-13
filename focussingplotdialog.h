#ifndef FOCUSSINGPLOTDIALOG_H
#define FOCUSSINGPLOTDIALOG_H

#include <vector>
#include <QDialog>
#include <QList>
#include <QPixmap>
#include <QVector>
#include <QString>

#include <qwt_plot_curve.h>
#include <qwt_symbol.h>

#include "modelfunction.h"

#include "ui_focussingCurvePlot.h"

#define FOCUSSINGPLOTDIALOG_NPOINTS 100

class FocusWidget; // just forward declaration

class FocussingPlotDialog : public QDialog
{
public:
    FocussingPlotDialog(const QString &title=QString::null, FocusWidget *parent = nullptr);

    void plot(QVector<double> &foc_value, QVector<double> &xfwhm, QVector<double> &yfwhm,
              size_t fit_order,
              QVector<double> &fit_coeffs);

private:
    Ui::FocussingPlotForm ui;

    QwtPlotCurve *xFWHMCurve, *yFWHMCurve;
    QwtPlotCurve *xFWHMFitCurve, *yFWHMFitCurve;

    FocusWidget *caller;

    QVector<double> focusValue, fwhmX, fwhmY;
    Polynomial focusRelation;
};

#endif // FOCUSSINGPLOTDIALOG_H
