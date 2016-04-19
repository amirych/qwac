#include "focuswidget.h"

#include <QDebug>


//#define FOCUSWIDGET_IMAGEPOINT_FMT "<b>X:</b> %1  <b>Y:</b> %2  <b>Value:</b> %3"
#define FOCUSWIDGET_IMAGEPOINT_FMT "<b>Pixel: </b> [%1; %2]  <b>  Value:</b> %3"

FocusWidget::FocusWidget(double start_val, double stop_val, double step_val, QWidget *parent): QMainWindow(parent),
    focusValueValidator(), focusImages(QStringList())
{
    ui.setupUi(this);

    imagePointLabel = new QLabel(QString(FOCUSWIDGET_IMAGEPOINT_FMT).arg(0.0,0,'f',1).arg(0.0,0,'f',1).arg(0.0,0,'f',1));
    this->statusBar()->addPermanentWidget(imagePointLabel);

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

