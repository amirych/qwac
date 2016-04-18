#include "focussingplotdialog.h"

#include "qwt_plot_curve.h"
#include <qwt_symbol.h>
#include <cmath>

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
    data_curve->setSymbol(QwtSymbol::Diamond);

    data_curve->setSamples(foc_value,fwhm);

    data_curve->attach(ui.focusCurvePlot);

    // fit curve (polynom)
    // coeffs: c0 + c1*x + c2*x^2 + ...

    QwtPlotCurve *fit_curve = new QwtPlotCurve;
    fit_curve->setStyle(QwtPlotCurve::Lines);

    data_curve->setSymbol(QwtSymbol::NoSymbol);

    QVector<double> fit_fwhm(fwhm.size(),fit_coeffs[0]);

    for ( size_t order = 1; order < fit_coeffs.size(); ++order ) {
        for ( size_t i = 0; i < foc_value.size(); ++i ) {
            fit_fwhm[i] += fit_coeffs[order]*std::pow(foc_value[i],order);
        }
    }
}
