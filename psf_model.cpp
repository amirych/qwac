#include "psf_model.h"

#include <limits>
#include <QDebug>

#define PSF_MODEL_FWHM_DEFAULT_MIN 0
#define PSF_MODEL_FWHM_DEFAULT_MAX std::numeric_limits<double>::infinity()

#define PSF_MODEL_CENTER_DEFAULT_MIN -std::numeric_limits<double>::infinity()
#define PSF_MODEL_CENTER_DEFAULT_MAX std::numeric_limits<double>::infinity()

// amplitude of PSF (maximal value)
#define PSF_MODEL_AMP_DEFAULT_MIN -std::numeric_limits<double>::infinity()
#define PSF_MODEL_AMP_DEFAULT_MAX std::numeric_limits<double>::infinity()

// Moffat's alpha-parameter
#define PSF_MODEL_ALPHA_DEFAULT_MIN 0.0
#define PSF_MODEL_ALPHA_DEFAULT_MAX std::numeric_limits<double>::infinity()

// model function rotation
#define PSF_MODEL_PHI_DEFAULT_MIN -std::numeric_limits<double>::infinity()
#define PSF_MODEL_PHI_DEFAULT_MAX std::numeric_limits<double>::infinity()


#define PSF_MODEL_FWHM_DEFAULT 3.0
#define PSF_MODEL_AMP_DEFAULT 0.0
#define PSF_MODEL_CENTER_DEFAULT 0.0
#define PSF_MODEL_ALPHA_DEFAULT 2.5
#define PSF_MODEL_PHI_DEFAULT 0.0


                /*  PSF model functions */

static void moffat(double *pars, double *func, int n_pars, int n_func, void* data);
static void gauss(double *pars, double *func, int n_pars, int n_func, void* data);

static void moffat2d(double *pars, double *func, int n_pars, int n_func, void* data);
static void gauss2d(double *pars, double *func, int n_pars, int n_func, void* data);



static QVector<double> emptyVector;


        /*  PSF_Model_params class realization  */

PSF_Model_Params::PSF_Model_Params(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub):
//    params(emptyVector), lowerBounds(emptyVector), upperBounds(emptyVector)
  params(pars), lowerBounds(lb), upperBounds(ub)
{
//    setParams(pars);
//    setLowerBounds(lb);
//    setUpperBounds(ub);
}


PSF_Model_Params::PSF_Model_Params(QVector<double> &pars):
    PSF_Model_Params(pars,emptyVector,emptyVector)
{
    setParams(pars);
}


PSF_Model_Params::PSF_Model_Params():
    PSF_Model_Params(emptyVector,emptyVector,emptyVector)
{

}


void PSF_Model_Params::setParams(QVector<double> &pars)
{
    params = pars;
}


void PSF_Model_Params::setContrains(QVector<double> &lb, QVector<double> &ub)
{
    setLowerBounds(lb);
    setUpperBounds(ub);
}

void PSF_Model_Params::setLowerBounds(QVector<double> &lb)
{
    qDebug() << "base Model lower bounds!";
    lowerBounds = lb;
}


void PSF_Model_Params::setUpperBounds(QVector<double> &ub)
{
    upperBounds = ub;
}


QVector<double> PSF_Model_Params::getParams() const
{
    return params;
}


QVector<double> PSF_Model_Params::getConstrains() const
{
    QVector<double> c = lowerBounds;
    c.resize(lowerBounds.size() + upperBounds.size());
    for ( int i = 0; i < upperBounds.size(); ++i) {
        c[lowerBounds.size() + i] = upperBounds[i];
    }

    return c;
}


QVector<double> PSF_Model_Params::getLowerBounds() const
{
    return lowerBounds;
}


QVector<double> PSF_Model_Params::getUpperBounds() const
{
    return upperBounds;
}



            /*  Gauss_Model_Params class realization  */

Gauss_Model_Params::Gauss_Model_Params(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub): PSF_Model_Params(pars,lb,ub)
{
    setParams(pars);
    setLowerBounds(lb);
    setUpperBounds(ub);
}


Gauss_Model_Params::Gauss_Model_Params(QVector<double> &pars): Gauss_Model_Params(pars,emptyVector,emptyVector)
{

}


Gauss_Model_Params::Gauss_Model_Params(): Gauss_Model_Params(emptyVector,emptyVector,emptyVector)
{

}

void Gauss_Model_Params::setParams(QVector<double> &pars)
{
    params = pars;

    if ( params.size() < PSF_MODEL_GAUSS_MIN_NPARS ) {
        int N = params.size();
        params.resize(PSF_MODEL_GAUSS_MIN_NPARS);
        int i = PSF_MODEL_GAUSS_MIN_NPARS;
        if ( N < i ) params[--i] = PSF_MODEL_FWHM_DEFAULT;
        if ( N < i ) params[--i] = PSF_MODEL_CENTER_DEFAULT;
        if ( N == 0 ) params[0] = PSF_MODEL_AMP_DEFAULT;
    }
}

void Gauss_Model_Params::setLowerBounds(QVector<double> &lb)
{
    qDebug() << "GAUSS lower bounds!";
    lowerBounds = lb;

    if ( lowerBounds.size() < PSF_MODEL_GAUSS_MIN_NPARS ) { // set default constrains
        int N = lowerBounds.size();
        lowerBounds.resize(PSF_MODEL_GAUSS_MIN_NPARS);
        int i = PSF_MODEL_GAUSS_MIN_NPARS;
        if ( N < i ) lowerBounds[--i] = PSF_MODEL_FWHM_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = PSF_MODEL_CENTER_DEFAULT_MIN;
        if ( N == 0 ) lowerBounds[0] = PSF_MODEL_AMP_DEFAULT_MIN;
    }
}


void Gauss_Model_Params::setUpperBounds(QVector<double> &ub)
{
    upperBounds = ub;
    if ( upperBounds.size() < PSF_MODEL_GAUSS_MIN_NPARS ) { // set default constrains
        int N = upperBounds.size();
        upperBounds.resize(PSF_MODEL_GAUSS_MIN_NPARS);
        int i = PSF_MODEL_GAUSS_MIN_NPARS;
        if ( N < i ) upperBounds[--i] = PSF_MODEL_FWHM_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = PSF_MODEL_CENTER_DEFAULT_MAX;
        if ( N == 0 ) upperBounds[0] = PSF_MODEL_AMP_DEFAULT_MAX;
    }
}

            /*  Gauss2D_Model_Params class realization  */

Gauss2D_Model_Params::Gauss2D_Model_Params(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub): PSF_Model_Params(pars,lb,ub)
{
    setParams(pars);
    setLowerBounds(lb);
    setUpperBounds(ub);
}


Gauss2D_Model_Params::Gauss2D_Model_Params(QVector<double> &pars): Gauss2D_Model_Params(pars,emptyVector,emptyVector)
{

}


Gauss2D_Model_Params::Gauss2D_Model_Params(): Gauss2D_Model_Params(emptyVector,emptyVector,emptyVector)
{

}

void Gauss2D_Model_Params::setParams(QVector<double> &pars)
{
    params = pars;

    if ( params.size() < PSF_MODEL_GAUSS2D_MIN_NPARS ) {
        int N = params.size();
        params.resize(PSF_MODEL_GAUSS2D_MIN_NPARS);
        int i = PSF_MODEL_GAUSS2D_MIN_NPARS;
        if ( N < i ) params[--i] = PSF_MODEL_PHI_DEFAULT;
        if ( N < i ) params[--i] = PSF_MODEL_FWHM_DEFAULT;
        if ( N < i ) params[--i] = PSF_MODEL_FWHM_DEFAULT;
        if ( N < i ) params[--i] = PSF_MODEL_CENTER_DEFAULT;
        if ( N < i ) params[--i] = PSF_MODEL_CENTER_DEFAULT;
        if ( N == 0 ) params[0] = PSF_MODEL_AMP_DEFAULT;
    }
}

void Gauss2D_Model_Params::setLowerBounds(QVector<double> &lb)
{
    qDebug() << "GAUSS2D lower bounds!";
    lowerBounds = lb;

    if ( lowerBounds.size() < PSF_MODEL_GAUSS2D_MIN_NPARS ) { // set default constrains
        int N = lowerBounds.size();
        lowerBounds.resize(PSF_MODEL_GAUSS2D_MIN_NPARS);
        int i = PSF_MODEL_GAUSS2D_MIN_NPARS;
        if ( N < i ) lowerBounds[--i] = PSF_MODEL_PHI_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = PSF_MODEL_FWHM_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = PSF_MODEL_FWHM_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = PSF_MODEL_CENTER_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = PSF_MODEL_CENTER_DEFAULT_MIN;
        if ( N == 0 ) lowerBounds[0] = PSF_MODEL_AMP_DEFAULT_MIN;
    }
}


void Gauss2D_Model_Params::setUpperBounds(QVector<double> &ub)
{
    upperBounds = ub;
    if ( upperBounds.size() < PSF_MODEL_GAUSS2D_MIN_NPARS ) { // set default constrains
        int N = upperBounds.size();
        upperBounds.resize(PSF_MODEL_GAUSS2D_MIN_NPARS);
        int i = PSF_MODEL_GAUSS2D_MIN_NPARS;
        if ( N < i ) upperBounds[--i] = PSF_MODEL_PHI_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = PSF_MODEL_FWHM_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = PSF_MODEL_FWHM_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = PSF_MODEL_CENTER_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = PSF_MODEL_CENTER_DEFAULT_MAX;
        if ( N == 0 ) upperBounds[0] = PSF_MODEL_AMP_DEFAULT_MAX;
    }
}



                    /*  Moffat_Model_Params class realization  */

Moffat_Model_Params::Moffat_Model_Params(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub): PSF_Model_Params(pars,lb,ub)
{
    setParams(pars);
    setLowerBounds(lb);
    setUpperBounds(ub);
}


Moffat_Model_Params::Moffat_Model_Params(QVector<double> &pars): Moffat_Model_Params(pars,emptyVector,emptyVector)
{

}


Moffat_Model_Params::Moffat_Model_Params(): Moffat_Model_Params(emptyVector,emptyVector,emptyVector)
{

}

void Moffat_Model_Params::setParams(QVector<double> &pars)
{
    params = pars;

    if ( params.size() < PSF_MODEL_MOFFAT_MIN_NPARS ) {
        int N = params.size();
        params.resize(PSF_MODEL_MOFFAT_MIN_NPARS);
        int i = PSF_MODEL_MOFFAT_MIN_NPARS;
        if ( N < i ) params[--i] = PSF_MODEL_ALPHA_DEFAULT;
        if ( N < i ) params[--i] = PSF_MODEL_FWHM_DEFAULT;
        if ( N < i ) params[--i] = PSF_MODEL_CENTER_DEFAULT;
        if ( N == 0 ) params[0] = PSF_MODEL_AMP_DEFAULT;
    }
}

void Moffat_Model_Params::setLowerBounds(QVector<double> &lb)
{
    qDebug() << "Moffat lower bounds!";
    lowerBounds = lb;

    if ( lowerBounds.size() < PSF_MODEL_MOFFAT_MIN_NPARS ) { // set default constrains
        int N = lowerBounds.size();
        lowerBounds.resize(PSF_MODEL_MOFFAT_MIN_NPARS);
        int i = PSF_MODEL_MOFFAT_MIN_NPARS;
        if ( N < i ) lowerBounds[--i] = PSF_MODEL_ALPHA_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = PSF_MODEL_FWHM_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = PSF_MODEL_CENTER_DEFAULT_MIN;
        if ( N == 0 ) lowerBounds[0] = PSF_MODEL_AMP_DEFAULT_MIN;
    }
}


void Moffat_Model_Params::setUpperBounds(QVector<double> &ub)
{
    upperBounds = ub;
    if ( upperBounds.size() < PSF_MODEL_MOFFAT_MIN_NPARS ) { // set default constrains
        int N = upperBounds.size();
        upperBounds.resize(PSF_MODEL_MOFFAT_MIN_NPARS);
        int i = PSF_MODEL_MOFFAT_MIN_NPARS;
        if ( N < i ) upperBounds[--i] = PSF_MODEL_ALPHA_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = PSF_MODEL_FWHM_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = PSF_MODEL_CENTER_DEFAULT_MAX;
        if ( N == 0 ) upperBounds[0] = PSF_MODEL_AMP_DEFAULT_MAX;
    }
}

            /*  Moffat2D_Model_Params class realization  */

Moffat2D_Model_Params::Moffat2D_Model_Params(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub): PSF_Model_Params(pars,lb,ub)
{
    setParams(pars);
    setLowerBounds(lb);
    setUpperBounds(ub);
}


Moffat2D_Model_Params::Moffat2D_Model_Params(QVector<double> &pars): Moffat2D_Model_Params(pars,emptyVector,emptyVector)
{

}


Moffat2D_Model_Params::Moffat2D_Model_Params(): Moffat2D_Model_Params(emptyVector,emptyVector,emptyVector)
{

}

void Moffat2D_Model_Params::setParams(QVector<double> &pars)
{
    params = pars;

    if ( params.size() < PSF_MODEL_MOFFAT2D_MIN_NPARS ) {
        int N = params.size();
        params.resize(PSF_MODEL_MOFFAT2D_MIN_NPARS);
        int i = PSF_MODEL_MOFFAT2D_MIN_NPARS;
        if ( N < i ) params[--i] = PSF_MODEL_PHI_DEFAULT;
        if ( N < i ) params[--i] = PSF_MODEL_ALPHA_DEFAULT;
        if ( N < i ) params[--i] = PSF_MODEL_FWHM_DEFAULT;
        if ( N < i ) params[--i] = PSF_MODEL_FWHM_DEFAULT;
        if ( N < i ) params[--i] = PSF_MODEL_CENTER_DEFAULT;
        if ( N < i ) params[--i] = PSF_MODEL_CENTER_DEFAULT;
        if ( N == 0 ) params[0] = PSF_MODEL_AMP_DEFAULT;
    }
}

void Moffat2D_Model_Params::setLowerBounds(QVector<double> &lb)
{
    qDebug() << "MoffatD lower bounds!";
    lowerBounds = lb;

    if ( lowerBounds.size() < PSF_MODEL_MOFFAT2D_MIN_NPARS ) { // set default constrains
        int N = lowerBounds.size();
        lowerBounds.resize(PSF_MODEL_MOFFAT2D_MIN_NPARS);
        int i = PSF_MODEL_MOFFAT2D_MIN_NPARS;
        if ( N < i ) lowerBounds[--i] = PSF_MODEL_PHI_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = PSF_MODEL_ALPHA_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = PSF_MODEL_FWHM_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = PSF_MODEL_FWHM_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = PSF_MODEL_CENTER_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = PSF_MODEL_CENTER_DEFAULT_MIN;
        if ( N == 0 ) lowerBounds[0] = PSF_MODEL_AMP_DEFAULT_MIN;
    }
}


void Moffat2D_Model_Params::setUpperBounds(QVector<double> &ub)
{
    upperBounds = ub;
    if ( upperBounds.size() < PSF_MODEL_MOFFAT2D_MIN_NPARS ) { // set default constrains
        int N = upperBounds.size();
        upperBounds.resize(PSF_MODEL_MOFFAT2D_MIN_NPARS);
        int i = PSF_MODEL_MOFFAT2D_MIN_NPARS;
        if ( N < i ) upperBounds[--i] = PSF_MODEL_PHI_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = PSF_MODEL_ALPHA_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = PSF_MODEL_FWHM_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = PSF_MODEL_FWHM_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = PSF_MODEL_CENTER_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = PSF_MODEL_CENTER_DEFAULT_MAX;
        if ( N == 0 ) upperBounds[0] = PSF_MODEL_AMP_DEFAULT_MAX;
    }
}


            /*  PSF_Model realization   */

PSF_Model::PSF_Model(QObject *parent) : QObject(parent)
{

}
