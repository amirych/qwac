#ifndef FOCUSSINGPLOTDIALOG_H
#define FOCUSSINGPLOTDIALOG_H

#include <QDialog>
#include <QList>
#include <QPixmap>
#include <QVector>
#include <QString>

#include "ui_focussingCurvePlot.h"


#define FOCUSSINGPLOTDIALOG_NPOINTS 100


class FocussingPlotDialog : public QDialog
{
public:
    FocussingPlotDialog(const QString &title=QString::null, QWidget *parent = nullptr);

    void plot(QVector<double> &foc_value, QVector<double> &fwhm,
              QVector<double> &fit_coeffs, QList<QPixmap> &images);

private:
    Ui::FocussingPlotForm ui;

};

#endif // FOCUSSINGPLOTDIALOG_H
