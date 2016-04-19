#include "focussingplotdialog.h"

#include "qwt_plot_curve.h"
#include <qwt_symbol.h>

FocussingPlotDialog::FocussingPlotDialog(const QString &title, QWidget *parent): QDialog(parent)
{
    ui.setupUi(this);

    ui.focusCurvePlot->setTitle(title);
    ui.focusCurvePlot->setAxisTitle(QwtPlot::yLeft,"FWHM");
    ui.focusCurvePlot->setAxisTitle(QwtPlot::xBottom,"Focus value");
}


void FocussingPlotDialog::plot(QVector<double> &foc_value, QVector<double> &fwhm,
                               QVector<double> &fit_coeffs, QList<QPixmap> &images)
{


    // measured point
    QwtPlotCurve *data_curve = new QwtPlotCurve;
    data_curve->setStyle(QwtPlotCurve::NoCurve);
    data_curve->setSymbol(new QwtSymbol(QwtSymbol::Diamond));

    data_curve->setSamples(foc_value,fwhm);

    data_curve->attach(ui.focusCurvePlot);

    // fit curve (polynom)
    // coeffs: c0 + c1*x + c2*x^2 + ...

    QwtPlotCurve *fit_curve = new QwtPlotCurve;
    fit_curve->setStyle(QwtPlotCurve::Lines);

    data_curve->setSymbol(new QwtSymbol(QwtSymbol::NoSymbol));

    QwtInterval foc_range = ui.focusCurvePlot->axisInterval(QwtPlot::xBottom);
    double foc_step = (foc_range.maxValue() - foc_range.minValue())/(FOCUSSINGPLOTDIALOG_NPOINTS-1);

    QVector<double> foc(FOCUSSINGPLOTDIALOG_NPOINTS);
    QVector<double> fit_fwhm(FOCUSSINGPLOTDIALOG_NPOINTS,fit_coeffs[0]);

    for ( size_t i = 0; i < FOCUSSINGPLOTDIALOG_NPOINTS; ++i ) {
        foc[i] = foc_range.minValue() + i*foc_step;
        double x = foc[i];
        for ( size_t order = 1; order < fit_coeffs.size(); ++order ) {
            fit_fwhm[i] += fit_coeffs[order]*x;
            x *= x;
        }
    }

    fit_curve->setSamples(foc,fit_fwhm);

    fit_curve->attach(ui.focusCurvePlot);
}
