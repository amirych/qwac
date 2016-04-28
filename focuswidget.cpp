#include "focuswidget.h"

#include <cmath>
#include <QDebug>
#include <QtConcurrent>

//#define FOCUSWIDGET_IMAGEPOINT_FMT "<b>X:</b> %1  <b>Y:</b> %2  <b>Value:</b> %3"
#define FOCUSWIDGET_IMAGEPOINT_FMT "<b>Pixel: </b> [%1; %2]  <b>  Value:</b> %3"
#define FOCUSWIDGET_MOVE_FOCUS_MSG_FMT "<b>Move focus to %1 position ...</b>"
#define FOCUSWIDGET_GET_IMAGE_MSG_FMT "<b>Get image: %1 ...</b>"
#define FOCUSWIDGET_OK_MSG "<b>  OK</b>"
#define FOCUSWIDGET_FAILURE_MSG "<b>Operation failed! (exit code: [%1])</b>"


FocusWidget::FocusWidget(double start_val, double stop_val, double step_val, QWidget *parent): QMainWindow(parent),
    currentFocusValue(0.0),
    focusValueValidator(), focusImages(QStringList()), focusPos(QVector<double>())
{
    ui.setupUi(this);

    imagePointLabel = new QLabel(QString(FOCUSWIDGET_IMAGEPOINT_FMT).arg(0.0,0,'f',1).arg(0.0,0,'f',1).arg(0.0,0,'f',1));
    this->statusBar()->addPermanentWidget(imagePointLabel);

    statusLabel = new QLabel("");
    this->statusBar()->addWidget(statusLabel);

//    this->setMouseTracking(true);
    ui.viewFrame->setFocusProxy(ui.view);

//    QVBoxLayout *view_layout = new QVBoxLayout(ui.fitsViewerFrame);
//    view_layout->setMargin(1);
//    view = new FitsViewWidget(ui.fitsViewerFrame);
//    view_layout->addWidget(view);
//    view->load("/home/timur/OLD_CYLON/WORK/s110808/S8950314.FTS");
//    view->showImage();

//    qDebug() << ui.view->getError();
//    qDebug() << ui.view->getZoom();

    setInitSetup(start_val, stop_val, step_val);

    ui.startLineEdit->setValidator(&focusValueValidator);
    ui.stopLineEdit->setValidator(&focusValueValidator);
    ui.stepLineEdit->setValidator(&focusValueValidator);

    expParamsDialog = new ExpParamsDialog(this);

    connect(ui.quitButton,SIGNAL(clicked(bool)),this,SLOT(about_quit()));
    connect(ui.runButton,SIGNAL(clicked(bool)),this,SLOT(run()));
    connect(ui.stopButton,SIGNAL(clicked(bool)),this,SLOT(stop()));
    connect(ui.focussingButton,SIGNAL(clicked(bool)),this,SLOT(focussing()));
    connect(ui.expPropsButton,SIGNAL(clicked(bool)),this,SLOT(setExpProps()));

    connect(ui.view,SIGNAL(imageIsShown(QString)),ui.filenameLineEdit,SLOT(setText(QString)));
    connect(ui.view,SIGNAL(imagePoint(QPointF,double)),this,SLOT(showImagePoint(QPointF,double)));

    focussingSequenceThread = new runSequence(this);

    connect(focussingSequenceThread,SIGNAL(started()),this, SLOT(sequenceIsStarting()));
    connect(focussingSequenceThread,SIGNAL(finished()),this,SLOT(sequenceIsFinished()));
    connect(focussingSequenceThread,SIGNAL(status(QString)),this,SLOT(setStatusMsg(QString)));

    // default root filename for image series
    QString filename = "foc";
    QStringList rates;
    rates << "Normal";
    expParamsDialog->init(filename,rates,0,1,1);
}

FocusWidget::FocusWidget(QWidget *parent): FocusWidget(0.0,0.0,0.0,parent)
{
}


FocusWidget::~FocusWidget()
{
    if ( focussingSequenceThread->isRunning() ) focussingSequenceThread->requestInterruption();
}

            /*   PUBLIC METHODS   */

void FocusWidget::setInitSetup(double start_val, double stop_val, double step_val)
{
    QString val;

    val.setNum(start_val);
    ui.startLineEdit->setText(val);

    val.setNum(stop_val);
    ui.stopLineEdit->setText(val);

    val.setNum(step_val);
    ui.stepLineEdit->setText(val);

    startFocusValue = start_val;
    stopFocusValue = stop_val;
    stepFocusValue = step_val;
}


void FocusWidget::setFocusValueRange(double min_val, double max_val, int decimals)
{
    focusValueValidator.setRange(min_val, max_val, decimals);
}


void FocusWidget::setExpInitSetup(QString &rootfilename, QStringList &rate, int rate_index, int xbin, int ybin)
{
    expParamsDialog->init(rootfilename,rate,rate_index,xbin,ybin);
}



            /*   PRIVATE SLOTS   */

void FocusWidget::showImagePoint(QPointF pos, double val)
{
    QString str = QString(FOCUSWIDGET_IMAGEPOINT_FMT).arg(pos.x(),0,'f',1).arg(pos.y(),0,'f',1).arg(val,0,'f',1);
    imagePointLabel->setText(str);
}


void FocusWidget::sequenceIsStarting()
{
    ui.focussingButton->setEnabled(false);
    ui.focussingButton->setToolTip("The focussing sequence is still in progress");
    ui.expPropsButton->setEnabled(false);
}


void FocusWidget::sequenceIsFinished()
{
    ui.focussingButton->setEnabled(true);
    ui.focussingButton->setToolTip("Select an object and press the button");
    ui.expPropsButton->setEnabled(true);
}


            /*   PUBLIC SLOTS   */

void FocusWidget::setStatusMsg(QString msg)
{
    statusLabel->setText(msg);
}


            /*   PROTECTED SLOTS   */

// base realization of reaction on "Quit" button (just close the widget)
void FocusWidget::about_quit()
{
    this->deleteLater();
}


void FocusWidget::run()
{
    startFocusValue = ui.startLineEdit->text().toDouble();
    stopFocusValue = ui.stopLineEdit->text().toDouble();
    stepFocusValue = ui.stepLineEdit->text().toDouble();

    double n_focus_pos = floor((stopFocusValue - startFocusValue)/stepFocusValue);

    focusPos.clear();
    focusImages.clear();

    int digit_num = floor(log10(n_focus_pos)+1.0);

    QString root_filename = expParamsDialog->getFilename();
    double exp_time = expParamsDialog->getExptime();

    for ( int i = 1; i <= n_focus_pos; ++i ) {
        focusPos << startFocusValue + (i-1)*stepFocusValue;
        focusImages << QString(root_filename + "%1.fits").arg(i,digit_num,10,QChar('0'));
    }

    focussingSequenceThread->initSequence(focusPos,focusImages,exp_time,nullptr);

    focussingSequenceThread->start();

//    QString status_str;
//    for ( int i = 0; i < focusPos.size(); ++i ) {
//        status_str = QString(FOCUSWIDGET_MOVE_FOCUS_MSG_FMT).arg(focusPos[i],0,'f');
//        statusLabel->setText(status_str);
//        int ret = moveFocus(focusPos[i]);
//        if ( ret ) {
//            status_str = QString(FOCUSWIDGET_FAILURE_MSG).arg(ret);
//            statusLabel->setText(status_str);
//            return;
//        } else {
//            statusLabel->setText(FOCUSWIDGET_OK_MSG);
//        }

//        status_str = QString(FOCUSWIDGET_GET_IMAGE_MSG_FMT).arg(focusImages[i]);
//        statusLabel->setText(status_str);
//        ret = getImage(focusImages[i],exp_time);
//        if ( ret ) {
//            status_str = QString(FOCUSWIDGET_FAILURE_MSG).arg(ret);
//            statusLabel->setText(status_str);
//            return;
//        } else {
//            statusLabel->setText(FOCUSWIDGET_OK_MSG);
//        }
//    }
}


void FocusWidget::stop()
{
    focussingSequenceThread->requestInterruption();
    ui.view->load("/home/timur/OLD_CYLON/WORK/s110808/S8950314.FTS");
    ui.view->showImage();
}


void FocusWidget::focussing()
{

}


void FocusWidget::setExpProps()
{
    expParamsDialog->exec();
}




runSequence::runSequence(FocusWidget *parent): QThread(parent),
    focusPos(QVector<double>()), focusImages(QStringList()), exp_time(0.0), expParams(nullptr)
{
    caller = parent;
}


void runSequence::initSequence(QVector<double> &focuspos, QStringList &images, double exptime, void *exp_pars)
{
    focusPos = focuspos;
    focusImages = images;
    exp_time = exptime;
    expParams = exp_pars;
}


void runSequence::run()
{
    QString status_str;
    for ( int i = 0; i < focusPos.size(); ++i ) {
        if ( isInterruptionRequested() ) return;

        status_str = QString(FOCUSWIDGET_MOVE_FOCUS_MSG_FMT).arg(focusPos[i],0,'f',1);
        emit status(status_str);
        int ret = caller->moveFocus(focusPos[i]);
        if ( ret ) {
            status_str = QString(FOCUSWIDGET_FAILURE_MSG).arg(ret);
            emit status(status_str);
            return;
        } else {
            status_str = FOCUSWIDGET_OK_MSG;
            emit status(status_str);
        }

        if ( isInterruptionRequested() ) return;

        status_str = QString(FOCUSWIDGET_GET_IMAGE_MSG_FMT).arg(focusImages[i]);
        emit status(status_str);
        ret = caller->getImage(focusImages[i],exp_time, expParams);
        if ( ret ) {
            status_str = QString(FOCUSWIDGET_FAILURE_MSG).arg(ret);
            emit status(status_str);
            return;
        } else {
            status_str = FOCUSWIDGET_OK_MSG;
            emit status(status_str);
        }
    }
}
