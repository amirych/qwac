#include "focussingplotdialog.h"

#include <algorithm>

FocussingPlotDialog::FocussingPlotDialog(const QString &title, QWidget *parent): QDialog(parent),
    xFWHMCurve(nullptr), yFWHMCurve(nullptr), xFWHMFitCurve(nullptr), yFWHMFitCurve(nullptr)
{
    ui.setupUi(this);

    ui.focusCurvePlotX->setTitle(title);
    ui.focusCurvePlotX->setAxisTitle(QwtPlot::yLeft,"FWHM along X-axis");
    ui.focusCurvePlotX->setAxisTitle(QwtPlot::xBottom,"Focus value");

    ui.focusCurvePlotY->setAxisTitle(QwtPlot::yLeft,"FWHM along Y-axis");
    ui.focusCurvePlotY->setAxisTitle(QwtPlot::xBottom,"Focus value");

    // measured point along X
    xFWHMCurve = new QwtPlotCurve;
    xFWHMCurve->setStyle(QwtPlotCurve::NoCurve);
    xFWHMCurve->setSymbol(new QwtSymbol(QwtSymbol::Diamond));

    // measured point along Y
    yFWHMCurve = new QwtPlotCurve;
    yFWHMCurve->setStyle(QwtPlotCurve::NoCurve);
    yFWHMCurve->setSymbol(new QwtSymbol(QwtSymbol::Diamond));

    // fit curves (polynom)
    // coeffs: c0 + c1*x + c2*x^2 + ...

    xFWHMFitCurve = new QwtPlotCurve;
    xFWHMFitCurve->setStyle(QwtPlotCurve::Lines);
    xFWHMFitCurve->setSymbol(new QwtSymbol(QwtSymbol::NoSymbol));

    yFWHMFitCurve = new QwtPlotCurve;
    yFWHMFitCurve->setStyle(QwtPlotCurve::Lines);
    yFWHMFitCurve->setSymbol(new QwtSymbol(QwtSymbol::NoSymbol));
}


void FocussingPlotDialog::plot(QVector<double> &foc_value, QVector<double> &xfwhm, QVector<double> &yfwhm,
                               QVector<double> &fit_coeffs, QList<QPixmap> &images)
{

    // clear possible previous plotted curves
    ui.focusCurvePlotX->detachItems(QwtPlotItem::Rtti_PlotCurve, false);

    xFWHMCurve->setSamples(foc_value,xfwhm);
    xFWHMCurve->attach(ui.focusCurvePlotX);


    yFWHMCurve->setSamples(foc_value,yfwhm);
    yFWHMCurve->attach(ui.focusCurvePlotY);


    auto foc_range = std::minmax_element(foc_value.constBegin(),foc_value.constEnd());

    double foc_step = (*foc_range.second - *foc_range.first)/(FOCUSSINGPLOTDIALOG_NPOINTS-1);

    int fit_order = fit_coeffs.size()/2; // fit_coeffs consists of polynom coefficients for X and Y
    QVector<double> foc(FOCUSSINGPLOTDIALOG_NPOINTS);
    QVector<double> fit_fwhmX(FOCUSSINGPLOTDIALOG_NPOINTS,fit_coeffs[0]);
    QVector<double> fit_fwhmY(FOCUSSINGPLOTDIALOG_NPOINTS,fit_coeffs[fit_order]);

    for ( size_t i = 0; i < FOCUSSINGPLOTDIALOG_NPOINTS; ++i ) {
        foc[i] = *foc_range.first + i*foc_step;
        double x = foc[i];
        for ( int order = 1; order < fit_order; ++order ) {
            fit_fwhmX[i] += fit_coeffs[order]*x;
            fit_fwhmY[i] += fit_coeffs[fit_order+order]*x;
            x *= x;
        }
    }

    xFWHMFitCurve->setSamples(foc,fit_fwhmX);
    yFWHMFitCurve->setSamples(foc,fit_fwhmY);

    xFWHMFitCurve->attach(ui.focusCurvePlotX);
    yFWHMFitCurve->attach(ui.focusCurvePlotY);

    ui.focusCurvePlotX->replot();
    ui.focusCurvePlotY->replot();

    // align axis scales
    double min_val, max_val;

    QwtInterval xfoc_range = ui.focusCurvePlotX->axisInterval(QwtPlot::xBottom);
    QwtInterval yfoc_range = ui.focusCurvePlotY->axisInterval(QwtPlot::xBottom);

    min_val = xfoc_range.minValue() <= yfoc_range.minValue() ? xfoc_range.minValue() : yfoc_range.minValue();
    max_val = xfoc_range.maxValue() >= yfoc_range.maxValue() ? xfoc_range.maxValue() : yfoc_range.maxValue();

    ui.focusCurvePlotX->setAxisScale(QwtPlot::xBottom,min_val,max_val);
    ui.focusCurvePlotY->setAxisScale(QwtPlot::xBottom,min_val,max_val);

    xfoc_range = ui.focusCurvePlotX->axisInterval(QwtPlot::yLeft);
    yfoc_range = ui.focusCurvePlotY->axisInterval(QwtPlot::yLeft);

    min_val = xfoc_range.minValue() <= yfoc_range.minValue() ? xfoc_range.minValue() : yfoc_range.minValue();
    max_val = xfoc_range.maxValue() >= yfoc_range.maxValue() ? xfoc_range.maxValue() : yfoc_range.maxValue();

    ui.focusCurvePlotX->setAxisScale(QwtPlot::yLeft,min_val,max_val);
    ui.focusCurvePlotY->setAxisScale(QwtPlot::yLeft,min_val,max_val);
}
