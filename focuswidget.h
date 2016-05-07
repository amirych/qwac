#ifndef FOCUSWIDGET_H
#define FOCUSWIDGET_H

#include "focuswidget_global.h"
#include "ui_focuswidget.h"

#include "FitsViewWidget.h"
#include "expparamsdialog.h"
#include "psf_model.h"

#include <QString>
#include <QDoubleValidator>
#include <QPointer>
#include <QStringList>
#include <QVector>
#include <QThread>

static QVector<double> empty_vector;


class runSequence; // just forward declaration

class runFitting: public QThread
{
    Q_OBJECT
public:
    enum FitError{OK, InvalidFilename, MemoryAllocationError, FitsError};
    runFitting(QWidget *parent);
    void initFitting(PSF_Model *psf_model, QStringList &foc_images, QRectF &fit_area);
    void run();
signals:
    void fittingParams(QVector<double> params);
    void fittingComplete();
    void error(runFitting::FitError err);
private:
    QStringList focusImages;
    QRectF fitArea;
    PSF_Model* psfModel;
};



//
//  Pure virtual class FocusWidget
//
class FOCUSWIDGETSHARED_EXPORT FocusWidget: public QMainWindow
{

    Q_OBJECT

public:
    FocusWidget(QWidget *parent = 0);
    FocusWidget(QString &root_failename, double start_val, double stop_val, double step_val, QWidget *parent = 0);

    ~FocusWidget();

    void setInitSetup(QString &root_failename, double start_val, double stop_val, double step_val);
    void setExpInitSetup(QString &rootfilename, QStringList &rate, int rate_index, int xbin, int ybin);
    void setFittingSetup(QVector<double> &init_pars, QVector<double> &lb = empty_vector, QVector<double> &ub = empty_vector);

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
    void setSelectedArea(QRectF area);
    void clearSelectedArea();
    void fittingComplete();
    void fittingParams(QVector<double> params);
    void fittingError(runFitting::FitError err);

private:
    Ui::FocusWidgetForm ui;

    ExpParamsDialog *expParamsDialog;

    QLabel *imagePointLabel;
    QLabel *statusLabel;

    double startFocusValue;
    double stopFocusValue;
    double stepFocusValue;
    double currentFocusValue;
//    QDoubleValidator focusValueValidator;

    QStringList focusImages;
    QVector<double> focusPos;

    PSF_Model *psfModel;
    QVector<double> fitParams;
    QVector<double> fitLowerBounds;
    QVector<double> fitUpperBounds;
    QVector<double> fitFWHM;
    QVector<double> fitFWHMCoeffs;

    friend class runSequence;

    runSequence *focussingSequenceThread;
    QRectF selectedArea;

    runFitting *runFittingThread;
};



class runSequence: public QThread
{
    Q_OBJECT
public:
    runSequence(FocusWidget* parent);
    void run();
    void initSequence(QVector<double> &focuspos, QStringList &images, double exptime, void *exp_pars);
    int sequenceLength() const; // return number of successfuly obtained images
signals:
    void status(QString msg);
    void imageIsReady(QString filename);
private:
    FocusWidget *caller;
    QVector<double> focusPos;
    QStringList focusImages;
    double exp_time;
    void *expParams;

    int n_images;
};


#endif // FOCUSWIDGET_H
