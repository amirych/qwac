#ifndef FOCUSWIDGET_H
#define FOCUSWIDGET_H

#include "focuswidget_global.h"
#include "ui_focuswidget.h"

#include "FitsViewWidget.h"
#include "expparamsdialog.h"

#include <QString>
#include <QDoubleValidator>
#include <QPointer>
#include <QStringList>

class FOCUSWIDGETSHARED_EXPORT FocusWidget: public QMainWindow
{

    Q_OBJECT

public:
    FocusWidget(QWidget *parent = 0);
    FocusWidget(double start_val, double stop_val, double step_val, QWidget *parent = 0);

    ~FocusWidget();

    void setInitSetup(double start_val, double stop_val, double step_val);
    void setExpInitSetup(QString &rootfilename, QStringList &rate, int rate_index, int xbin, int ybin);

    void setFocusValueRange(double min_val, double max_val, int decimals = 1);

protected slots:
    virtual void about_quit();
    virtual void run();
    virtual void stop();
    virtual void focussing();
    virtual void setExpProps();

//    void resizeEvent(QResizeEvent *event);

//protected:
//    virtual int getImage(QString filename, double exp_time, void *exp_params = nullptr) = 0;
//    virtual int moveFocus(double focus_value) = 0;


private slots:
    void showImagePoint(QPointF pos, double val);

private:
    Ui::FocusWidgetForm ui;

    ExpParamsDialog *expParamsDialog;

    QLabel *imagePointLabel;

    double startFocusValue;
    double stopFocusValue;
    double stepFocusValue;
    QDoubleValidator focusValueValidator;

    QStringList focusImages;
};

#endif // FOCUSWIDGET_H
