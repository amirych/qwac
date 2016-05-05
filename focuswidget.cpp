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


#define FOCUSWIDGET_FWHM_FITTING_DEGREE 3 // parabola

// parabola function definition (synopsis is according to levmar package requirements)
// ext_data is interpreted as pointer to double array with independant variable (must be of length not lesser than n_func)
// pars is vector of coefficients of polynom in form [a0,a1,a2] and y = a0 + a1*x + a2*x^2
static void parabola_func(double* pars, double* func, int n_params, int n_func, void *ext_data)
{
    if ( ext_data == nullptr || pars = nullptr || func == nullptr ) return;
    if ( n_params <= 0 || n_func <= 0 ) return;

    double* x = (double*) ext_data;

    for ( int i = 0; i < n_func; ++i ) {
        func[i] = pars[0] + pars[1]*x[i] + pars[2]*x[i]*x[i];
    }
}



FocusWidget::FocusWidget(double start_val, double stop_val, double step_val, QWidget *parent): QMainWindow(parent),
    currentFocusValue(0.0),
    focusValueValidator(), focusImages(QStringList()), focusPos(QVector<double>()),
    psfModel(nullptr),
    fitParams(empty_vector), fitLowerBounds(empty_vector), fitUpperBounds(empty_vector),
    fitFWHM(empty_vector), fitFWHMCoeffs(QVector<double>(FOCUSWIDGET_FWHM_FITTING_DEGREE)),
    selectedArea(QRectF(0,0,0,0))
{
    ui.setupUi(this);

    imagePointLabel = new QLabel(QString(FOCUSWIDGET_IMAGEPOINT_FMT).arg(0.0,0,'f',1).arg(0.0,0,'f',1).arg(0.0,0,'f',1));
    this->statusBar()->addPermanentWidget(imagePointLabel);

    statusLabel = new QLabel("");
    this->statusBar()->addWidget(statusLabel);

//    this->setMouseTracking(true);
    ui.viewFrame->setFocusProxy(ui.view);

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

    connect(focussingSequenceThread,SIGNAL(imageIsReady(QString)),ui.view,SLOT(load(QString)));

    ui.focussingButton->setEnabled(false); // focussing button is not available still
    connect(ui.view,SIGNAL(regionWasDeselected()),this,SLOT(clearSelectedArea()));

    runFittingThread = new runFitting(this);
    connect(runFittingThread,SIGNAL(fittingComplete()),this,SLOT(fittingComplete()));
    connect(runFittingThread,SIGNAL(fittingParams(QVector<double>)),this,SLOT(fittingParams(QVector<double>)));
    connect(runFittingThread,SIGNAL(error(runFitting::FitError)),this,SLOT(fittingError(runFitting::FitError)));

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

    delete psfModel;
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
    // fit parabola to focus-FWHM relation

    QVector<double> lb(FOCUSWIDGET_FWHM_FITTING_DEGREE);
    QVector<double> ub(FOCUSWIDGET_FWHM_FITTING_DEGREE);

    QVector xFWHM, yFWHM;

    // set constrains
    for ( int i = 0; i < FOCUSWIDGET_FWHM_FITTING_DEGREE; ++i ) {
        lb[i] = -std::numeric_limits<double>::infinity();
        ub[i] = std::numeric_limits<double>::infinity();
    }

    // polynom coefficient of x^2 (for parabola) term must be greater than 0
    lb[FOCUSWIDGET_FWHM_FITTING_DEGREE-1] = 0.0;

    // set initial coefficients
    fitFWHMCoeffs[FOCUSWIDGET_FWHM_FITTING_DEGREE-1] = 1.0;
    for ( int i = 0; i < FOCUSWIDGET_FWHM_FITTING_DEGREE-2; ++i ) fitFWHMCoeffs[i] = 0.0;

    double *work_space = nullptr;

    int ret;

    try {
        work_space = new double[LM_BC_DIF_WORKSZ(3,focusPos.size())];
        for ( int i = 0; i < focusPos.size(); ++i ) {
            xFWHM.append(fitFWHM[i*2]);
            yFWHM.append(fitFWHM[i*2+1]);
        }
    } catch (std::bad_alloc &ex) {
        setStatusMsg("Memory allocation error occured!!!");
        return;
    }


    ret = dlevmar_bc_dif(parabola_func,fitFWHMCoeffs.data(),xFWHM.data(),FOCUSWIDGET_FWHM_FITTING_DEGREE,xFWHM.size()),
                         lb.data(),ub.data(),NULL,100,NULL,NULL,work_space,NULL,(void*)focusPos.data());

    delete[] work_space;

    if ( ret < 0 ) {
        setStatusMsg("Cannot fit focusPos - FWHM relation!!!");
        return;
    }

    // plot the results

}


void FocusWidget::fittingParams(QVector<double> params)
{
    fitFWHM.append(params[3]); // FWHM along X-axis
    fitFWHM.append(params[4]); // FWHM along Y-axis
}


void FocusWidget::fittingError(runFitting::FitError err)
{
    ui.focussingButton->setEnabled(true);

    QString str = QString("Fitting proccess failed! Error code: %1").arg((int)err);

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
    startFocusValue = ui.startLineEdit->text().toDouble();
    stopFocusValue = ui.stopLineEdit->text().toDouble();
    stepFocusValue = ui.stepLineEdit->text().toDouble();

    if ( stepFocusValue <= 0 ) {
        if ( startFocusValue >= stopFocusValue ) return;
    } else {
        if ( startFocusValue >= stopFocusValue ) return;
    }

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
}


void FocusWidget::stop()
{
    focussingSequenceThread->requestInterruption();
    ui.view->load("/home/timur/OLD_CYLON/WORK/s110808/S8950314.FTS");
    ui.view->showImage();
}


void FocusWidget::focussing()
{
    if ( selectedArea.isNull() ) return;

    ui.focussingButton->setEnabled(false);

    setStatusMsg("PSF model fitting ...");

    fitFWHM.clear();

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
        }
    }
}



runFitting::runFitting(QWidget *parent): QThread(parent),
  focusImages(QStringList()), fitArea(QRectF(0,0,0,0)),
  psfModel(nullptr)
{

}


void runFitting::initFitting(PSF_Model *psf_model, QStringList &foc_images, QRectF &fit_area)
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

    long xmin,ymin,xmax,ymax;

    xmin = (long) round(fitArea.x());
    ymin = (long) round(fitArea.y());
    xmax = (long) round(fitArea.x()+fitArea.width()-1);
    ymax = (long) round(fitArea.y()+fitArea.height()-1);

    QVector<double> x,y;
    for ( long xx = xmin; xx <= xmax; ++xx ) {
        for ( long yy = ymin; yy <= ymax; ++yy ) {
            x.append(xx);
            y.append(yy);
        }
    }

    psfModel->setArgument(x,y);

    for ( int i = 0; i < focusImages.length(); ++i ) {
        QString filename = focusImages[i].trimmed();
        if ( filename.isEmpty() ) {
            emit error(runFitting::InvalidFilename);
            return;
        }

        // special filename form for sub-image reading
        filename = QString(filename+"[%1:%2, %3:%4]").arg(xmin).arg(xmax).arg(ymin).arg(ymax);
        char* fname = filename.toLocal8Bit().data();

        try {
            fits_open_image(&FITS_fptr, fname, READONLY, &fits_status);
            if ( fits_status ) throw fits_status;

            fits_read_imghdr(FITS_fptr, maxdim, NULL, &bitpix, &naxis, naxes, NULL, NULL, NULL, &fits_status);
            if ( fits_status ) throw fits_status;

            for ( int i = 0; i < maxdim; ++i ) {
                nelem *= naxes[i];
            }

            image = std::unique_ptr<double[]>(new double[nelem]);
            double *buffer = image.get();

            fits_read_img(FITS_fptr, TDOUBLE, 1, nelem, NULL, (void*) buffer, NULL, &fits_status);
            if ( fits_status ) throw fits_status;

            fits_close_file(FITS_fptr, &fits_status);
            if ( fits_status ) throw fits_status;

            std::vector<double> v = std::vector<double>(buffer,buffer+nelem);
            QVector<double> data = QVector<double>::fromStdVector(v);

            psfModel->fitData(data);

            emit fittingParams(psfModel->getParams());

        } catch (std::bad_alloc &ex) {
            image = nullptr;
            fits_close_file(FITS_fptr, &fits_status);
            emit error(runFitting::MemoryAllocationError);
            return;
        } catch (int err) {
            image = nullptr;
            fits_close_file(FITS_fptr, &fits_status);
            emit error(runFitting::FitsError);
            return;
        }

    }

    emit fittingComplete();
}
