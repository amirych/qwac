#ifndef FOCUSSINGPLOTDIALOG_H
#define FOCUSSINGPLOTDIALOG_H

#include <QDialog>
#include <QList>
#include <QPixmap>
#include <QVector>
#include <QString>

#include <qwt_plot_curve.h>
#include <qwt_symbol.h>

#include "ui_focussingCurvePlot.h"

#define FOCUSSINGPLOTDIALOG_NPOINTS 100

class FocusWidget; // jst forward declaration

class FocussingPlotDialog : public QDialog
{
public:
    FocussingPlotDialog(const QString &title=QString::null, FocusWidget *parent = nullptr);

    void plot(QVector<double> &foc_value, QVector<double> &xfwhm, QVector<double> &yfwhm,
              QVector<double> &fit_coeffs);

private:
    Ui::FocussingPlotForm ui;

    QwtPlotCurve *xFWHMCurve, *yFWHMCurve;
    QwtPlotCurve *xFWHMFitCurve, *yFWHMFitCurve;

    FocusWidget *caller;

};

#endif // FOCUSSINGPLOTDIALOG_H
