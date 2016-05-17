#ifndef FOCUSSINGPLOTDIALOG_H
#define FOCUSSINGPLOTDIALOG_H

#include <vector>
#include <QDialog>
#include <QList>
#include <QPixmap>
#include <QVector>
#include <QString>
#include <QHBoxLayout>

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
    ~FocussingPlotDialog();

    void plot(QVector<double> &foc_value, QVector<double> &xfwhm, QVector<double> &yfwhm,
              size_t fit_order, QStringList &images, QRectF &view_area);

private:
    Ui::FocussingPlotForm ui;

    QwtPlotCurve *xFWHMCurve, *yFWHMCurve;
    QwtPlotCurve *xFWHMFitCurve, *yFWHMFitCurve;

    FocusWidget *caller;
    QHBoxLayout *images_layout;

    QVector<double> focusValue, fwhmX, fwhmY;
    Polynomial focusRelation;
};

#endif // FOCUSSINGPLOTDIALOG_H
