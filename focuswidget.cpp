#include "focuswidget.h"
//#include "psf_model.h"

#include <cmath>
#include <QDebug>
#include <QtConcurrent>
#include <levmar.h>
#include <fitsio.h>

//#define FOCUSWIDGET_IMAGEPOINT_FMT "<b>X:</b> %1  <b>Y:</b> %2  <b>Value:</b> %3"
#define FOCUSWIDGET_IMAGEPOINT_FMT "<b>Pixel: </b> [%1; %2]  <b>  Value:</b> %3"
#define FOCUSWIDGET_MOVE_FOCUS_MSG_FMT "<b>Move focus to %1 position ...</b>"
#define FOCUSWIDGET_GET_IMAGE_MSG_FMT "<b>Get image: %1 ...</b>"
#define FOCUSWIDGET_OK_MSG "<b>  OK</b>"
#define FOCUSWIDGET_FAILURE_MSG "<b>Operation failed! (exit code: [%1])</b>"

#define FOCUSWIDGET_FWHM_FITTING_ORDER 3 // parabola

static QString default_rootfilename = "foc";

// parabola function definition (synopsis is according to levmar package requirements)
// ext_data is interpreted as pointer to double array with independant variable (must be of length not lesser than n_func)
// pars is vector of coefficients of polynom in form [a0,a1,a2] and y = a0 + a1*x + a2*x^2
static void parabola_func(double* pars, double* func, int n_params, int n_func, void *ext_data)
{
    if ( ext_data == nullptr || pars == nullptr || func == nullptr ) return;
    if ( n_params <= 0 || n_func <= 0 ) return;

    double* x = (double*) ext_data;

    for ( int i = 0; i < n_func; ++i ) {
        func[i] = pars[0] + pars[1]*x[i] + pars[2]*x[i]*x[i];
    }
}


// generate PSF fitting model arguments
//static QVector<double> generate_arguments(QRectF &area, QVector<double> &x, QVector<double> &y) {
//    double xmin,ymin,xmax,ymax;

//    xmin = round(area.x());
//    ymin = round(area.y());
//    xmax = round(area.x()+area.width()-1);
//    ymax = round(area.y()+area.height()-1);

//    x.clear();
//    y.clear();
//    for ( double yy = ymin; yy <= ymax; ++yy ) {
//        for ( double xx = xmin; xx <= xmax; ++xx ) {
//            x.append(xx);
//            y.append(yy);
//        }
//    }

//    return QVector<double>({xmin,xmax,ymin,ymax});
//}


FocusWidget::FocusWidget(QString root_filename, double start_val, double stop_val, double step_val, QWidget *parent): QMainWindow(parent),
    currentFocusValue(0.0),
    focusImages(QStringList()), focusPos(QVector<double>()),
    psfModel(nullptr),
    fitParams(empty_vector), fitLowerBounds(empty_vector), fitUpperBounds(empty_vector),
    fitFWHM(empty_vector), fitFWHMCoeffs(QVector<double>(2*FOCUSWIDGET_FWHM_FITTING_ORDER)),
    psfModelX(std::vector<double>()), psfModelY(std::vector<double>()),
    selectedArea(QRectF(0,0,0,0))
{
    ui.setupUi(this);

    imagePointLabel = new QLabel(QString(FOCUSWIDGET_IMAGEPOINT_FMT).arg(0.0,0,'f',1).arg(0.0,0,'f',1).arg(0.0,0,'f',1));
    this->statusBar()->addPermanentWidget(imagePointLabel);

    statusLabel = new QLabel("");
    this->statusBar()->addWidget(statusLabel);

//    this->setMouseTracking(true);
    ui.viewFrame->setFocusProxy(ui.view);

    setInitSetup(root_filename, start_val, stop_val, step_val);

//    ui.startLineEdit->setValidator(&focusValueValidator);
//    ui.stopLineEdit->setValidator(&focusValueValidator);
//    ui.stepLineEdit->setValidator(&focusValueValidator);

    expParamsDialog = new ExpParamsDialog(this);

    connect(ui.rootnameLineEdit,SIGNAL(textChanged(QString)),expParamsDialog,SLOT(setFilename(QString)));

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

    connect(focussingSequenceThread,SIGNAL(imageIsReady(QString)),ui.view,SLOT(load(QString)));

    ui.focussingButton->setEnabled(false); // focussing button is not available still
    connect(ui.view,SIGNAL(regionWasDeselected()),this,SLOT(clearSelectedArea()));

    runFittingThread = new runFitting(this);
    connect(runFittingThread,SIGNAL(fittingComplete()),this,SLOT(fittingComplete()));
//    connect(runFittingThread,SIGNAL(fittingParams(QVector<double>)),this,SLOT(fittingParams(QVector<double>)));
    connect(runFittingThread,SIGNAL(error(int)),this,SLOT(fittingError(int)));

    connect(runFittingThread, &runFitting::started,
          [=]{ ui.focussingButton->setEnabled(false); ui.runButton->setEnabled(false);});
//    connect(runFittingThread, &runFitting::finished,
//          [=]{ ui.focussingButton->setEnabled(true); ui.runButton->setEnabled(true); });
//    connect(runFittingThread, &runFitting::fittingComplete,
//          [=]{ ui.focussingButton->setEnabled(true); ui.runButton->setEnabled(true); });


    psfModelBackgroundDegree[0] = 0;
    psfModelBackgroundDegree[1] = 0;

    // defaults for image series
    QStringList rates;
    rates << "Normal";

    expParamsDialog->init(root_filename,rates,0,1,1);

    plotDialog = new FocussingPlotDialog("",this);
}

FocusWidget::FocusWidget(QWidget *parent): FocusWidget(default_rootfilename, 0.0,0.0,0.0,parent)
{
}


FocusWidget::~FocusWidget()
{
    if ( focussingSequenceThread->isRunning() ) focussingSequenceThread->requestInterruption();

    delete psfModel;

//    delete focussingSequenceThread;
//    delete runFittingThread;
//    delete plotDialog;
}

            /*   PUBLIC METHODS   */

void FocusWidget::setInitSetup(QString &root_filename, double start_val, double stop_val, double step_val)
{
    ui.rootnameLineEdit->setText(root_filename);
    ui.startFocusSpinBox->setValue(start_val);

    ui.stopFocusSpinBox->setValue(stop_val);

    ui.stepFocusSpinBox->setValue(step_val);

    startFocusValue = start_val;
    stopFocusValue = stop_val;
    stepFocusValue = step_val;
}


void FocusWidget::setFocusValueRange(double min_val, double max_val, int decimals)
{
//    focusValueValidator.setRange(min_val, max_val, decimals);

    ui.startFocusSpinBox->setRange(min_val, max_val);
    ui.startFocusSpinBox->setDecimals(decimals);

    ui.stopFocusSpinBox->setRange(min_val, max_val);
    ui.stopFocusSpinBox->setDecimals(decimals);

    ui.stepFocusSpinBox->setRange(min_val, max_val);
    ui.stepFocusSpinBox->setDecimals(decimals);
}


void FocusWidget::setExpInitSetup(QString &rootfilename, QStringList &rate, int rate_index, int xbin, int ybin)
{
    expParamsDialog->init(rootfilename,rate,rate_index,xbin,ybin);
}


void FocusWidget::setFittingSetup(QVector<double> &init_pars, QVector<double> &lb, QVector<double> &ub)
{
    fitParams = init_pars;
    fitLowerBounds = lb;
    fitUpperBounds = ub;
}


            /*   PRIVATE SLOTS   */

void FocusWidget::showImagePoint(QPointF pos, double val)
{
    QString str = QString(FOCUSWIDGET_IMAGEPOINT_FMT).arg(pos.x(),0,'f',1).arg(pos.y(),0,'f',1).arg(val,0,'f',1);
    imagePointLabel->setText(str);
}


void FocusWidget::sequenceIsStarting()
{
//    ui.focussingButton->setEnabled(false);
    ui.focussingButton->setToolTip("The focussing sequence is still in progress");
    ui.expPropsButton->setEnabled(false);
    disconnect(ui.view,SIGNAL(regionWasSelected(QRectF)),this,SLOT(setSelectedArea(QRectF)));
}


void FocusWidget::sequenceIsFinished()
{
//    ui.focussingButton->setEnabled(true);
    ui.focussingButton->setToolTip("Select an object and press the button");
    ui.expPropsButton->setEnabled(true);

    connect(ui.view,SIGNAL(regionWasSelected(QRectF)),this,SLOT(setSelectedArea(QRectF)));

    // show focussing image with possible near-optimal focus

    int N_images = focussingSequenceThread->sequenceLength();
    if ( N_images < focusPos.size() ) {
        if ( !N_images ) { // no one image
            setStatusMsg("No one image was obtained!!!");
            return;
        }

        for (int i = N_images; i < focusPos.length(); ++i ) focusImages.removeLast();
        focusPos.resize(N_images);
    }

    int N = focusPos.size()/2; // assume near-optimal focus is on image at middle of focussing images sequence

    ui.view->load(focusImages[N]);
//    qDebug() << focusImages[N] << "   (" << ui.view->getError() << ")";

    setStatusMsg("Select an object and push 'Focussing' button");
}


void FocusWidget::setSelectedArea(QRectF area)
{
    selectedArea = area;
    ui.focussingButton->setEnabled(true);
}

void FocusWidget::clearSelectedArea()
{
    selectedArea.setWidth(0.0);
    selectedArea.setHeight(0.0);
    ui.focussingButton->setEnabled(false);
}


void FocusWidget::fittingComplete()
{
    ui.runButton->setEnabled(true);

    setStatusMsg("Compute the best focus and plotting the results ...");

//    // fit parabola to focus-FWHM relation

//    QVector<double> lb(FOCUSWIDGET_FWHM_FITTING_ORDER);
//    QVector<double> ub(FOCUSWIDGET_FWHM_FITTING_ORDER);

//    QVector<double> xFWHM, yFWHM;

//    // set constrains
//    for ( int i = 0; i < FOCUSWIDGET_FWHM_FITTING_ORDER; ++i ) {
//        lb[i] = -std::numeric_limits<double>::infinity();
//        ub[i] = std::numeric_limits<double>::infinity();
//    }

//    // polynom coefficient of x^2 (for parabola) term must be greater than 0
//    lb[FOCUSWIDGET_FWHM_FITTING_ORDER-1] = 0.0;

//    // set initial coefficients
//    for ( int i = 0; i < 2*FOCUSWIDGET_FWHM_FITTING_ORDER; ++i ) fitFWHMCoeffs[i] = 0.0;
//    fitFWHMCoeffs[FOCUSWIDGET_FWHM_FITTING_ORDER-1] = 1.0; // along X-axis
//    fitFWHMCoeffs[2*FOCUSWIDGET_FWHM_FITTING_ORDER-1] = 1.0; // along Y-axis

//    double *work_space = nullptr;

//    int ret;

//    try {
//        work_space = new double[LM_BC_DIF_WORKSZ(FOCUSWIDGET_FWHM_FITTING_ORDER,focusPos.size())];
//        for ( int i = 0; i < focusPos.size(); ++i ) {
//            xFWHM.append(fitFWHM[i*2]);
//            yFWHM.append(fitFWHM[i*2+1]);
//        }
//    } catch (std::bad_alloc &ex) {
//        setStatusMsg("Memory allocation error occured!!!");
//        return;
//    }


//    // fit relation for FWHM along X-axis
//    ret = dlevmar_bc_dif(parabola_func,fitFWHMCoeffs.data(),xFWHM.data(),FOCUSWIDGET_FWHM_FITTING_ORDER,xFWHM.size(),
//                         lb.data(),ub.data(),NULL,100,NULL,NULL,work_space,NULL,(void*)focusPos.data());


//    if ( ret < 0 ) {
//        setStatusMsg("Cannot fit 'focusPos - FWHM' relation!!!");
//        delete[] work_space;
//        return;
//    }

//    // fit relation for FWHM along Y-axis
//    ret = dlevmar_bc_dif(parabola_func,fitFWHMCoeffs.data()+FOCUSWIDGET_FWHM_FITTING_ORDER,yFWHM.data(),FOCUSWIDGET_FWHM_FITTING_ORDER,yFWHM.size(),
//                         lb.data(),ub.data(),NULL,100,NULL,NULL,work_space,NULL,(void*)focusPos.data());


//    if ( ret < 0 ) {
//        setStatusMsg("Cannot fit 'focusPos - FWHM' relation!!!");
//        delete[] work_space;
//        return;
//    }

//    delete[] work_space;


    // plot the results

    QVector<std::vector<double>> pars;
    runFittingThread->getFitParams(pars);

    QVector<double> xFWHM, yFWHM;
    for ( int i = 0; i < focusPos.size(); ++i ) {
        xFWHM.append(pars[i][3]);
        yFWHM.append(pars[i][4]);
    }

    plotDialog->plot(focusPos,xFWHM,yFWHM,FOCUSWIDGET_FWHM_FITTING_ORDER,fitFWHMCoeffs);

    plotDialog->exec();
}


void FocusWidget::fittingParams(QVector<double> params)
{
    fitFWHM.append(params[3]); // FWHM along X-axis
    fitFWHM.append(params[4]); // FWHM along Y-axis
    qDebug() << "fitted params: " << params;
}


void FocusWidget::fittingError(int err)
{
    ui.focussingButton->setEnabled(true);

    QString str = QString("Fitting proccess failed! Error code: %1").arg(err);

    setStatusMsg(str);
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
    startFocusValue = ui.startFocusSpinBox->value();
    stopFocusValue = ui.stopFocusSpinBox->value();
    stepFocusValue = ui.stepFocusSpinBox->value();

    if ( stepFocusValue <= 0 ) {
        if ( startFocusValue >= stopFocusValue ) return;
    } else {
        if ( startFocusValue >= stopFocusValue ) return;
    }

    double n_focus_pos = floor((stopFocusValue - startFocusValue)/stepFocusValue);

    focusPos.clear();
    focusImages.clear();

    int digit_num = floor(log10(n_focus_pos)+2.0);

    QString root_filename = expParamsDialog->getFilename();
    double exp_time = expParamsDialog->getExptime();

    for ( int i = 0; i <= n_focus_pos; ++i ) {
        focusPos << startFocusValue + i*stepFocusValue;
        focusImages << QString(root_filename + "%1.fits").arg((i+1),digit_num,10,QChar('0'));
    }

    focussingSequenceThread->initSequence(focusPos,focusImages,exp_time,nullptr);

    qDebug() << "START SEQUENCE: " << focusPos;
    qDebug() << "START SEQUENCE: " << focusImages;
    qDebug() << "START SEQUENCE: " << focussingSequenceThread->isFinished();

    focussingSequenceThread->start();
}


void FocusWidget::stop()
{
    focussingSequenceThread->requestInterruption();
//    ui.view->load("/home/timur/foc04.fits");
//    ui.view->showImage();
}


void FocusWidget::focussing()
{
    if ( selectedArea.isNull() ) return;

    setStatusMsg("PSF model fitting ...");

    fitFWHM.clear();

//    QVector<double> x,y;

//    generate_arguments(selectedArea,x,y);

//    size_t back_deg[2] = {0,0}; // just constant as background
    std::vector<double> pp(3);
//    std::vector<double> pp(5);
    pp[0] = 0.0;
    pp[1] = selectedArea.x() + selectedArea.width()/2.0;
    pp[2] =  selectedArea.y() + selectedArea.height()/2.0;
//    pp[3] = 20;
//    pp[4] = 20;

    if ( ui.psfModelComboBox->currentText() == "Moffat2D" ) {
        psfModel = new Moffat2DModelFunction(pp,psfModelBackgroundDegree);
    } else if ( ui.psfModelComboBox->currentText() == "Gauss2D" ) {
        psfModel = new Gauss2DModelFunction(pp,psfModelBackgroundDegree);
    } else {
            setStatusMsg("Unknown PSF model!!!");
            return;
    }

//    psfModelX = x.toStdVector();
//    psfModelY = y.toStdVector();

//    psfModel->setArgument(psfModelX,psfModelY);

    double xmin = round(selectedArea.x());
    double ymin = round(selectedArea.y());
    double xmax = round(selectedArea.x() + selectedArea.width());
    double ymax = round(selectedArea.y() + selectedArea.height());

    psfModel->setArgument(xmin,xmax,ymin,ymax);

//    pp = psfModel->getParams();
//    pp[pp.size()-1] = 3022.0;
//    psfModel->setParams(pp,psfModelBackgroundDegree);

//    qDebug() << "ARG LEN: " << psfModelX.size();
    qDebug() << "INIT PARS: " << QVector<double>::fromStdVector(pp);
//    psfModel->setParams(pp,(void*)back_deg);

    // set initial center of PSF model as the center of the selected area
//    if ( fitParams.length() < 3) { // user does not setup initial parameters or it is invalid number of parameters
//        fitParams = psfModel->getParams();
//    }
//    fitParams[1] = selectedArea.x() + selectedArea.width()/2.0;
//    fitParams[2] = selectedArea.y() + selectedArea.height()/2.0;
//    psfModel->setParams(fitParams);


    runFittingThread->initFitting(psfModel,focusImages,selectedArea);

    runFittingThread->start();
}


void FocusWidget::setExpProps()
{
    expParamsDialog->exec();
}




runSequence::runSequence(FocusWidget *parent): QThread(parent),
    focusPos(QVector<double>()), focusImages(QStringList()), exp_time(0.0), expParams(nullptr), n_images(0)
{
    caller = parent;
}


void runSequence::initSequence(QVector<double> &focuspos, QStringList &images, double exptime, void *exp_pars)
{
    focusPos = focuspos;
    focusImages = images;
    exp_time = exptime;
    expParams = exp_pars;

    n_images = 0;
}


int runSequence::sequenceLength() const
{
    return n_images;
}

void runSequence::run()
{
    qDebug() << "RUN METHOD!";
    QString status_str;
    for ( n_images = 0; n_images < focusPos.size(); ++n_images ) {
        if ( isInterruptionRequested() ) { // still no read images, so -1 for n_image
            --n_images;
            return;
        }

        status_str = QString(FOCUSWIDGET_MOVE_FOCUS_MSG_FMT).arg(focusPos[n_images],0,'f',1);
        emit status(status_str);
        int ret = caller->moveFocus(focusPos[n_images]);
        if ( ret ) {
            status_str = QString(FOCUSWIDGET_FAILURE_MSG).arg(ret);
            emit status(status_str);
            --n_images; // still no read images, so -1 for n_image
            return;
        } else {
            status_str = FOCUSWIDGET_OK_MSG;
            emit status(status_str);
        }

        if ( isInterruptionRequested() ) { // still no read images, so -1 for n_image
            --n_images;
            return;
        }

        status_str = QString(FOCUSWIDGET_GET_IMAGE_MSG_FMT).arg(focusImages[n_images]);
        emit status(status_str);
        ret = caller->getImage(focusImages[n_images],exp_time, expParams);
        if ( ret ) {
            status_str = QString(FOCUSWIDGET_FAILURE_MSG).arg(ret);
            emit status(status_str);
            --n_images; // still no read images, so -1 for n_image
            return;
        } else {
            status_str = FOCUSWIDGET_OK_MSG;
            emit status(status_str);
            emit imageIsReady(focusImages[n_images]);
        }
    }
}



runFitting::runFitting(QWidget *parent): QThread(parent),
  focusImages(QStringList()), fitArea(QRectF(0,0,0,0)),
  psfModel(nullptr), fitParams(QVector<std::vector<double>>())
{

}


//void runFitting::initFitting(PSF_Model *psf_model, QStringList &foc_images, QRectF &fit_area)
void runFitting::initFitting(ModelFunction2D *psf_model, QStringList &foc_images, QRectF &fit_area)
{
    focusImages = foc_images;
    fitArea = fit_area;
    psfModel = psf_model;
}


void runFitting::run()
{
    if ( psfModel == nullptr ) return;
    if ( focusImages.length() == 0 ) return;
    if ( fitArea.isEmpty() ) return;

    fitsfile *FITS_fptr;
    int fits_status = 0;

    // THE ONLY 2D-images is support now!!!
    int maxdim = 2;
    long naxes[maxdim];
    int naxis, bitpix;
    LONGLONG nelem = 1;

    std::unique_ptr<double[]> image;


//    QVector<double> x,y;

//    QVector<double> borders = generate_arguments(fitArea,x,y);

    double xmin = round(fitArea.x());
    double ymin = round(fitArea.y());
    double xmax = round(fitArea.x() + fitArea.width());
    double ymax = round(fitArea.y() + fitArea.height());

//    psfModel->setArgument(x,y);

    for ( int i = 0; i < focusImages.length(); ++i ) {
        QString filename = focusImages[i].trimmed();
        if ( filename.isEmpty() ) {
            emit error((int) 104); // CFITSIO error
            return;
        }

        // special filename form for sub-image reading
        filename = QString(filename+"[%1:%2, %3:%4]").arg(xmin,0,'f',0).arg(xmax,0,'f',0).arg(ymin,0,'f',0).arg(ymax,0,'f',0);
//        filename = QString(filename+"[%1:%2, %3:%4]").arg(borders[0],0,'f',0).arg(borders[1],0,'f',0).arg(borders[2],0,'f',0).arg(borders[3],0,'f',0);
//        char* fname = filename.toLocal8Bit().data();

        char fname[FLEN_FILENAME];
        strncpy(fname,filename.toLocal8Bit().data(),filename.size()+1);

        qDebug() << "PSF fitting: " << filename;

        try {
            fits_open_image(&FITS_fptr, fname, READONLY, &fits_status);
            if ( fits_status ) throw fits_status;

            fits_read_imghdr(FITS_fptr, maxdim, NULL, &bitpix, &naxis, naxes, NULL, NULL, NULL, &fits_status);
            if ( fits_status ) throw fits_status;

            nelem = 1;
            for ( int i = 0; i < maxdim; ++i ) {
                nelem *= naxes[i];
            }

            image = std::unique_ptr<double[]>(new double[nelem]);
            double *buffer = image.get();

            fits_read_img(FITS_fptr, TDOUBLE, 1, nelem, NULL, (void*) buffer, NULL, &fits_status);
            if ( fits_status ) throw fits_status;

            fits_close_file(FITS_fptr, &fits_status);
            if ( fits_status ) throw fits_status;

            fits_status = 0;

//            std::vector<double> v = std::vector<double>(buffer,buffer+nelem);
//            QVector<double> data = QVector<double>::fromStdVector(v);

//            std::vector<double> data = std::vector<double>(buffer,buffer+nelem);

//            psfModel->fitData(data);

            psfModel->fitData(buffer,nelem);

            fitParams.append(psfModel->getParams());

//            std::vector<double> lb,ub;
//            psfModel->getConstrains(&lb,&ub);
//            qDebug() << "lower bounds: " << QVector<double>::fromStdVector(lb);
//            qDebug() << "upper bounds: " << QVector<double>::fromStdVector(ub);

//            qDebug() << "DATA Nelem: " << nelem;
//            qDebug() << "FIT PARS: " << QVector<double>::fromStdVector(psfModel->getParams());
//            qDebug() << "FIT STATUS: " << psfModel->getFitStatus();

//            emit fittingParams(psfModel->getParams());

        } catch (std::bad_alloc &ex) {
            throw (int) -1;
//            image = nullptr;
//            fits_close_file(FITS_fptr, &fits_status);
//            emit error(runFitting::MemoryAllocationError);
//            return;
        } catch (int err) {
            image = nullptr;
            fits_close_file(FITS_fptr, &fits_status);
//            qDebug() << "PSF fitting err: " << err;
//            emit error(runFitting::FitsError);
            emit error(err);
            return;
        }

    }

    emit fittingComplete();
}


void runFitting::getFitParams(QVector<std::vector<double> > &pars)
{
    pars.clear();
    pars = fitParams;
}
