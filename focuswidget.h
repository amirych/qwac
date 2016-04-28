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
#include <QVector>
#include <QThread>


class runSequence; // just forward declaration

//
//  Pure virtual class
//
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

public slots:
    void setStatusMsg(QString msg);

protected slots:
    virtual void about_quit();
    virtual void run();
    virtual void stop();
    virtual void focussing();
    virtual void setExpProps();

//    void resizeEvent(QResizeEvent *event);

protected:
//    virtual int getImage(QString filename, double exp_time, void *exp_params = nullptr) = 0;
//    virtual int moveFocus(double focus_value) = 0;
    virtual int getImage(QString filename, double exp_time, void *exp_params = nullptr) {QThread::sleep(3); return 0;}
    virtual int moveFocus(double focus_value) {QThread::sleep(3); return 0;}


private slots:
    void showImagePoint(QPointF pos, double val);
    void sequenceIsStarting();
    void sequenceIsFinished();

private:
    Ui::FocusWidgetForm ui;

    ExpParamsDialog *expParamsDialog;

    QLabel *imagePointLabel;
    QLabel *statusLabel;

    double startFocusValue;
    double stopFocusValue;
    double stepFocusValue;
    double currentFocusValue;
    QDoubleValidator focusValueValidator;

    QStringList focusImages;
    QVector<double> focusPos;

    friend class runSequence;

    runSequence *focussingSequenceThread;
};



class runSequence: public QThread
{
    Q_OBJECT
public:
    runSequence(FocusWidget* parent);
    void run();
    void initSequence(QVector<double> &focuspos, QStringList &images, double exptime, void *exp_pars);
signals:
    void status(QString msg);
private:
    FocusWidget *caller;
    QVector<double> focusPos;
    QStringList focusImages;
    double exp_time;
    void *expParams;
};

#endif // FOCUSWIDGET_H
