#include "focuswidget.h"

#include <cmath>
#include <QDebug>


//#define FOCUSWIDGET_IMAGEPOINT_FMT "<b>X:</b> %1  <b>Y:</b> %2  <b>Value:</b> %3"
#define FOCUSWIDGET_IMAGEPOINT_FMT "<b>Pixel: </b> [%1; %2]  <b>  Value:</b> %3"
#define FOCUSWIDGET_MOVE_FOCUS_MSG_FMT "<b>Move focus to [%1] position ...</b>"
#define FOCUSWIDGET_GET_IMAGE_MSG_FMT "<b>Get image: [%1] ...</b>"
#define FOCUSWIDGET_OK_MSG "<b>OK</b>"


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
        focusImages << QString(root_filename + "[%1]").arg(i,digit_num,10,QChar('0'));
    }

    for ( int i = 0; i < focusPos.size(); ++i ) {
        statusLabel->setText(QString(FOCUSWIDGET_MOVE_FOCUS_MSG_FMT).arg(focusPos[i]),0,'f');
        int ret = moveFocus(focusPos[i]);
        if ( !ret ) {

        } else {

        }
        ret = getImage(focusImages[i],exp_time);
    }
}


void FocusWidget::stop()
{
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

