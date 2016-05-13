#include "focussingplotdialog.h"

#include <algorithm>

#include "focuswidget.h"

FocussingPlotDialog::FocussingPlotDialog(const QString &title, FocusWidget *parent): QDialog((QWidget*)parent),
    xFWHMCurve(nullptr), yFWHMCurve(nullptr), xFWHMFitCurve(nullptr), yFWHMFitCurve(nullptr), caller(parent),
    fwhmX(QVector<double>()), fwhmY(QVector<double>()), focusRelation(Polynomial())
{
    ui.setupUi(this);


    ui.focusCurvePlotX->setTitle(title);
    ui.focusCurvePlotX->setAxisTitle(QwtPlot::yLeft,"FWHM along X-axis");
    ui.focusCurvePlotX->setAxisTitle(QwtPlot::xBottom,"Focus value");

    ui.focusCurvePlotY->setAxisTitle(QwtPlot::yLeft,"FWHM along Y-axis");
    ui.focusCurvePlotY->setAxisTitle(QwtPlot::xBottom,"Focus value");

    QPen sym_pen = QPen(QColor("black"));
    QBrush sym_brush = QBrush(Qt::SolidPattern);
//    QSize sym_size = QSize(5,5);
    QwtSymbol *curve_sym = new QwtSymbol(QwtSymbol::Rect);
    curve_sym->setBrush(sym_brush);
    curve_sym->setPen(sym_pen);
    curve_sym->setSize(7);

    // measured point along X
    xFWHMCurve = new QwtPlotCurve;
    xFWHMCurve->setStyle(QwtPlotCurve::NoCurve);
    xFWHMCurve->setSymbol(curve_sym);


    // measured point along Y
    yFWHMCurve = new QwtPlotCurve;
    yFWHMCurve->setStyle(QwtPlotCurve::NoCurve);
    yFWHMCurve->setSymbol(curve_sym);

    // fit curves (polynom)
    // coeffs: c0 + c1*x + c2*x^2 + ...

    QColor col("red");
    xFWHMFitCurve = new QwtPlotCurve;
    xFWHMFitCurve->setStyle(QwtPlotCurve::Lines);
    xFWHMFitCurve->setSymbol(new QwtSymbol(QwtSymbol::NoSymbol));
    xFWHMFitCurve->setPen(col,1.0);

    yFWHMFitCurve = new QwtPlotCurve;
    yFWHMFitCurve->setStyle(QwtPlotCurve::Lines);
    yFWHMFitCurve->setSymbol(new QwtSymbol(QwtSymbol::NoSymbol));
    yFWHMFitCurve->setPen(col,1.0);
}


void FocussingPlotDialog::plot(QVector<double> &foc_value, QVector<double> &xfwhm, QVector<double> &yfwhm, size_t fit_order,
                               QVector<double> &fit_coeffs)
{

    // clear possible previous plotted curves
    ui.focusCurvePlotX->detachItems(QwtPlotItem::Rtti_PlotCurve, false);

    xFWHMCurve->setSamples(foc_value,xfwhm);
    xFWHMCurve->attach(ui.focusCurvePlotX);


    yFWHMCurve->setSamples(foc_value,yfwhm);
    yFWHMCurve->attach(ui.focusCurvePlotY);

    auto foc_range = std::minmax_element(foc_value.constBegin(),foc_value.constEnd());

    qDebug() << "PARABOLA FITTING: ";
//    qDebug() << "focVal: " << foc_value;
//    qDebug() << "xfwhm: " << xfwhm;
//    qDebug() << "yfwhm: " << yfwhm;

    // fit "focus - FWHM" relations

    QVector<double> fit_coeffsX, fit_coeffsY;
    std::vector<double> coeffs(fit_order);
    std::vector<double> lb(fit_order),ub(fit_order);
    coeffs[0] = (*foc_range.second + *foc_range.first)/2.0;

    focusRelation.setParams(coeffs);
    focusRelation.getConstrains(&lb,&ub);
    lb[fit_order-1] = 0.0; // set lower constrains
    focusRelation.setConstrains(lb,ub);

    std::vector<double> arg = foc_value.toStdVector();
    focusRelation.setArgument(arg);

    arg = xfwhm.toStdVector();
    focusRelation.fitData(arg);
    fit_coeffsX = QVector<double>::fromStdVector(focusRelation.getParams());
//    qDebug() << "PARABOLA X: " << fit_coeffsX;

    arg = yfwhm.toStdVector();
    focusRelation.fitData(arg);
    fit_coeffsY = QVector<double>::fromStdVector(focusRelation.getParams());
//    qDebug() << "PARABOLA Y: " << fit_coeffsY;

    double foc_step = (*foc_range.second - *foc_range.first)/(FOCUSSINGPLOTDIALOG_NPOINTS-1);

//    int fit_order = fit_coeffs.size()/2; // fit_coeffs consists of polynom coefficients for X and Y
    QVector<double> foc(FOCUSSINGPLOTDIALOG_NPOINTS);
    QVector<double> fit_fwhmX(FOCUSSINGPLOTDIALOG_NPOINTS,fit_coeffsX[0]);
    QVector<double> fit_fwhmY(FOCUSSINGPLOTDIALOG_NPOINTS,fit_coeffsY[0]);

    for ( size_t i = 0; i < FOCUSSINGPLOTDIALOG_NPOINTS; ++i ) {
        foc[i] = *foc_range.first + i*foc_step;
        double x = foc[i];
        for ( int order = 1; order < fit_order; ++order ) {
            fit_fwhmX[i] += fit_coeffsX[order]*x;
            fit_fwhmY[i] += fit_coeffsY[order]*x;
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

    ui.focusCurvePlotX->replot();
    ui.focusCurvePlotY->replot();

    // compute the best focus value
    double best_focusX, best_focusY, best_focus;
    if ( fit_order == 3 ) {
        best_focusX = -fit_coeffsX[1]/fit_coeffsX[2]/2.0;
        best_focusY = -fit_coeffsY[1]/fit_coeffsY[2]/2.0;

        best_focus = ( best_focusX + best_focusY)/2.0; // just mean
    }

    ui.bestFocusSpinBox->setValue(best_focus);

    connect(ui.setBestFocusButton,&QPushButton::clicked,[=]{caller->moveFocus(best_focus);});
}
