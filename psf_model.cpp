#include "psf_model.h"

#include <limits>
#include <QDebug>

// number of model parameters without possible background polynom parameters
#define PSF_MODEL_GAUSS_MIN_NPARS 3      // [Amplitude, center, fwhm]
#define PSF_MODEL_MOFFAT_MIN_NPARS 4     // [Amplitude, center, fwhm, alpha]
#define PSF_MODEL_GAUSS2D_MIN_NPARS 6    // [Amplitude, centerX, centerY, fwhmX, fwhmY, phi]
#define PSF_MODEL_MOFFAT2D_MIN_NPARS 7   // [Amplitude, centerX, centerY, fwhmX, fwhmY, alpha, phi]


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

// polynom function coefficient default constrains
#define PSF_MODEL_POLY_DEFAULT_MIN -std::numeric_limits<double>::infinity()
#define PSF_MODEL_POLY_DEFAULT_MAX std::numeric_limits<double>::infinity()

// default constrains
#define PSF_MODEL_BOUND_DEFAULT_MIN -std::numeric_limits<double>::infinity()
#define PSF_MODEL_BOUND_DEFAULT_MAX std::numeric_limits<double>::infinity()

#define PSF_MODEL_FWHM_DEFAULT 3.0
#define PSF_MODEL_AMP_DEFAULT 0.0
#define PSF_MODEL_CENTER_DEFAULT 0.0
#define PSF_MODEL_ALPHA_DEFAULT 2.5
#define PSF_MODEL_PHI_DEFAULT 0.0


                /*  Polynom functions  */

//
//  p(x) = coeffs[0] + coeffs[1]*x + coeffs[2]*x^2 + ... + coeffs[n_coeffs-1]*x^(n_coeffs-1)
//
static void poly(double *x, double *coeffs, int n_x, int n_coeffs, double *val){
    double arg;

    for ( int i = 0; i < n_x; ++i ) {
        val[i] = coeffs[0];
        arg = x[i];
        for ( int ord = 1; ord < n_coeffs; ++ord, arg *= arg) {
            val[i] += coeffs[ord]*arg;
        }
    }
}

//
// p(x,y) = coefss[0] + coeffs[1]*y + coeffs[2]*y^2 + ... + coeffs[yd]*y^yd +
//        + coeffs[yd+1]*x + coeffs[yd+2]*x*y + coeffs[yd+3]*x*y^2 + ... + coeffs[2*(yd+1)-1]*x*y^yd +
//        + coeffs[2*(yd+1)]*x^2 + coeffs[2*(yd+1)+1]*x^2*y + coeffs[2*(yd+1)+2]*x^2*y^2 + ... + coeffs[3*(yd+1)-1]*x^2*y^yd +
//        ...
//        + coeffs[xd*(yd+1)]*x^xd + coeffs[xd*(yd+1)+1]*x^xd*y + coeffs[xd*(yd+1)+2]*x^xd*y^2 + ... + coeffs[(xd+1)*(y+1)-1]*x^xd*y^yd
//
//  , where xd and yd are polynom degree along X and Y
//
static void poly2d(double *x, double *y, double *coeffs, int n, int x_degree, int y_degree, double *val){
    double argx,argy;
    int idx = 1;

    for ( int i = 0; i < n; ++i ) {
        val[i] = coeffs[0];
        argx = 1.0;
        for ( int ord_x = 0; ord_x <= x_degree; ++ord_x, argx *= x[i] ) {
            argy = 1.0;
            for ( int ord_y = 0; ord_y <= y_degree; ++ord_y, argy *= y[i] ) {
                val[i] = coeffs[idx++]*argx*argy;
            }

        }
    }
}

                /*  PSF model functions */

//
//  Gaussian:
//    g(x) = a0*exp(-z^2/2.0), where z = (x-a1)/sigma and sigma = a2/2.0/sqrt(2.0*ln(2.0))
//    a0 - amplitude
//    a1 - center
//    a2 - FWHM
//
static int gauss(QVector<double> &pars, QVector<double> &func, extraData_t *data){
    if ( pars.size() == 0 ) return 1;
    if ( pars.size() < PSF_MODEL_GAUSS_MIN_NPARS) return 1;
    if ( data == nullptr ) return 1;

    if ( data->x.size() == 0 ) return 1;
    if ( data->x.size() != func.size() ) func.resize(data->x.size());

    int n_poly = pars.size() - PSF_MODEL_GAUSS_MIN_NPARS;

    func.fill(0.0);

    // compute background polynom
    if ( n_poly ) {
        poly(data->x.data(),pars.data()+PSF_MODEL_GAUSS_MIN_NPARS,data->x.size(),n_poly,func.data());
    }

    // compute gaussian
    double sigma = pars[2]/2.35482004503; // FWHM to gaussian sigma

    for ( int i = 0; i < func.size(); ++i ) {
        double z2 = (data->x[i] - pars[1])/sigma;
        z2 *= z2;
        func[i] += pars[0]*exp(-z2/2.0);
    }

    return 0;
}


//
//  Moffat function:
//    M(x) = a0*(1+z^2)^(-a4), where z = (x-a1)/s and s = a2/2/sqrt(2^(1/a4)-1)
//    a0 - amplitude
//    a1 - center
//    a2 - FWHM
//    a3 - moffat's parameter
//
static int moffat(QVector<double> &pars, QVector<double> &func, extraData_t *data){
    if ( pars.size() == 0 ) return 1;
    if ( pars.size() < PSF_MODEL_MOFFAT_MIN_NPARS) return 1;
    if ( data == nullptr ) return 1;

    if ( data->x.size() == 0 ) return 1;
    if ( data->x.size() != func.size() ) func.resize(data->x.size());

    int n_poly = pars.size() - PSF_MODEL_MOFFAT_MIN_NPARS;

    func.fill(0.0);

    // compute background polynom
    if ( n_poly ) {
        poly(data->x.data(),pars.data()+PSF_MODEL_MOFFAT_MIN_NPARS,data->x.size(),n_poly,func.data());
    }

    // compute Moffat function
    for ( int i = 0; i < func.size(); ++i ) {
        double term = pars[2]/2.0/sqrt(pow(2.0,1.0/pars[3])-1.0);
        double z2 = (data->x[i] - pars[1])/term;
        z2 *= z2;

        func[i] += pars[0]*pow(1.0+z2,-pars[3]);
    }

    return 0;
}

/*
    Bivariate Gaussian:

        g(x,y) = A*exp(-k*r^2), where

            k = 4*ln(2)
            r^2 = ((x-xc)/fwhmX)^2 + ((y-yc)/fwhmY)^2

             The present routine realization allows to compute Gaussian function
             with its symmetry axes that are not parallel to coordinate axes.
             Following this, in the first the routine computes transformation to new "tilted"
             coordinate system:

               x' = (x-xc)*cos(phi) - (y-yc)*sin(phi)
               y' = (y-yc)*cos(phi) + (x-xc)*sin(phi), where

               phi - rotation angle (clock-wise from X-axis)

     pars = [A, xc, yc, fwhmX, fwhmY, phi, ...], where ... denotes possible additive 2D-polynom coefficients
            phi is in radians

*/
static int gauss2d(QVector<double> &pars, QVector<double> &func, extraData_t *data){
    if ( pars.size() == 0 ) return 1;
    if ( pars.size() < PSF_MODEL_GAUSS2D_MIN_NPARS) return 1;
    if ( data == nullptr ) return 1;

    if ( data->x.size() == 0 ) return 1;
    if ( data->y.size() == 0 ) return 1;
    if ( data->x.size() != data->y.size() ) return 1;


    int n_poly;
    int *degree;

    if ( data->extra_params == nullptr ) n_poly = 0;
    else {
        n_poly = pars.size() - PSF_MODEL_GAUSS2D_MIN_NPARS;
        degree = (int*)data->extra_params; // should be at least 2-element int-type array!!!
        if ( n_poly < (degree[0]+1)*(degree[1]+1) ) return 1;
    }

    if ( data->x.size() != func.size() ) func.resize(data->x.size());
    func.fill(0.0);

    if ( n_poly ) { // compute 2D-polynom
        poly2d(data->x.data(),data->y.data(),pars.data()+ PSF_MODEL_GAUSS2D_MIN_NPARS,data->x.size(),
               degree[0],degree[1],func.data());
    }

    double k = 4.0*log(2.0);

    for ( int i = 0; i < func.size(); ++i ) {
        // first, rotate cordinate system
        double t = data->x[i] - pars[1];
        double y = data->y[i] - pars[2];

        double x = t*cos(pars[5]) - y*sin(pars[5]);
        y = y*cos(pars[5]) + t*sin(pars[5]);

        double r1 = x/pars[3];
        double r2 = y/pars[4];
        double z2 = r1*r1 + r2*r2;

        func[i] += pars[0]*exp(-k*z2);
    }

    return 0;
}


/*  The Moffat's bivariate function:

              f(x,y) = A/(1+k*r(x,y)^2)^alpha, where

               A - amplitude of the function
               k = 4*(2^(1/alpha)-1) - norm. parameter:
               r(x,y) = r^2 = ((x-xc)/fwhmX)^2+((y-yc)/fwhmY)^2 - radius vector
                 xc and yc - x and y peak positions
                 fwhmX and fwhmY - full width at half maximum along X and Y axes, respectively
               alpha - slope of Moffat's function

             The present routine realization allows to compute Moffat's function
             with its symmetry axes that are not parallel to coordinate axes.
             Following this, in the first the routine computes transformation to new "tilted"
             coordinate system:

               x' = (x-xc)*cos(phi) - (y-yc)*sin(phi)
               y' = (y-yc)*cos(phi) + (x-xc)*sin(phi), where

               phi - rotation angle (clock-wise from X-axis)

     pars = [A, xc, yc, fwhmX, fwhmY, alpha, phi, ...], where ... denotes possible additive 2D-polynom coefficients
            phi is in radians

*/
static int moffat2d(QVector<double> &pars, QVector<double> &func, extraData_t *data){
    if ( pars.size() == 0 ) return 1;
    if ( pars.size() < PSF_MODEL_MOFFAT2D_MIN_NPARS) return 1;
    if ( data == nullptr ) return 1;

    if ( data->x.size() == 0 ) return 1;
    if ( data->y.size() == 0 ) return 1;
    if ( data->x.size() != data->y.size() ) return 1;


    int n_poly;
    int *degree;

    if ( data->extra_params == nullptr ) n_poly = 0;
    else {
        n_poly = pars.size() - PSF_MODEL_MOFFAT2D_MIN_NPARS;
        degree = (int*)data->extra_params; // should be at least 2-element int-type array!!!
        if ( n_poly < (degree[0]+1)*(degree[1]+1) ) return 1;
    }

    if ( data->x.size() != func.size() ) func.resize(data->x.size());
    func.fill(0.0);

    if ( n_poly ) { // compute 2D-polynom
        poly2d(data->x.data(),data->y.data(),pars.data()+ PSF_MODEL_MOFFAT2D_MIN_NPARS,data->x.size(),
               degree[0],degree[1],func.data());
    }

    double k = 2.0*(pow(2.0,1.0/pars[5])-1.0);

    for ( int i = 0; i < func.size(); ++i ) {
        // first, rotate cordinate system
        double t = data->x[i] - pars[1];
        double y = data->y[i] - pars[2];

        double x = t*cos(pars[6]) - y*sin(pars[6]);
        y = y*cos(pars[6]) + t*sin(pars[6]);

        double r1 = x/pars[3];
        double r2 = y/pars[4];
        double z2 = r1*r1 + r2*r2;

        func[i] += pars[0]*pow(1.0+k*z2,-pars[5]);
    }

    return 0;
}


static QVector<double> emptyVector;


        /*  PSF_Model class realization  */


PSF_Model::PSF_Model(psfModelFunc_t func, QVector<double> &pars, extraData_t *extra_data):
    modelFunc(func), modelFuncExtraData(extra_data), params(pars), modelFuncValue(emptyVector)
{
    if ( extra_data == nullptr ) {
        modelFuncExtraData = new extraData_t;
        modelFuncExtraData->extra_params = nullptr;
    }
}


PSF_Model::PSF_Model(psfModelFunc_t func):
    PSF_Model(func,emptyVector,nullptr)
{

}


PSF_Model::~PSF_Model()
{
    delete modelFuncExtraData;
}


void PSF_Model::setParams(QVector<double> &pars, extraData_t *extra_data)
{
    params = pars;
    modelFuncExtraData = extra_data;
}


void PSF_Model::setContrains(QVector<double> &lb, QVector<double> &ub)
{
    setLowerBounds(lb);
    setUpperBounds(ub);
}

void PSF_Model::setLowerBounds(QVector<double> &lb)
{
    qDebug() << "base Model lower bounds!";
    lowerBounds = lb;
}


void PSF_Model::setUpperBounds(QVector<double> &ub)
{
    upperBounds = ub;
}


QVector<double> PSF_Model::getParams() const
{
    return params;
}


QVector<double> PSF_Model::getConstrains() const
{
    QVector<double> c = lowerBounds;
    c.resize(lowerBounds.size() + upperBounds.size());
    for ( int i = 0; i < upperBounds.size(); ++i) {
        c[lowerBounds.size() + i] = upperBounds[i];
    }

    return c;
}


QVector<double> PSF_Model::getLowerBounds() const
{
    return lowerBounds;
}


QVector<double> PSF_Model::getUpperBounds() const
{
    return upperBounds;
}


void PSF_Model::setArgument(QVector<double> &x)
{
    modelFuncExtraData->x = x;
}

void PSF_Model::setArgument(QVector<double> &x, QVector<double> &y)
{
    modelFuncExtraData->x = x;
    modelFuncExtraData->y = y;
}


QVector<double> PSF_Model::operator ()()
{
    compute();
    return modelFuncValue;
}


QVector<double> PSF_Model::operator ()(QVector<double> &x)
{
    modelFuncExtraData->x = x;

    compute();

    return modelFuncValue;
}


QVector<double> PSF_Model::operator ()(QVector<double> &x, QVector<double> &y)
{
    modelFuncExtraData->x = x;
    modelFuncExtraData->y = y;

    compute();
    return modelFuncValue;
}


void PSF_Model::ensureConstrains()
{
    int n = lowerBounds.size();
    if (  n < params.size() ) {
        lowerBounds.resize(params.size());
        for (int i = n; i < params.size(); ++i ) lowerBounds[i] = PSF_MODEL_BOUND_DEFAULT_MIN;
    }

    n = upperBounds.size();
    if (  n < params.size() ) {
        upperBounds.resize(params.size());
        for (int i = n; i < params.size(); ++i ) upperBounds[i] = PSF_MODEL_BOUND_DEFAULT_MAX;
    }
}


void PSF_Model::objetive_function(double *pars, double *func, int n_pars, int n_func, void *data)
{
    compute();
}


void PSF_Model::compute()
{
    int ret_code = modelFunc(params,modelFuncValue,modelFuncExtraData);
    if ( ret_code ) {

    }
}


            /*  Gauss_Model class realization  */

Gauss_Model::Gauss_Model(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub):
    PSF_Model(gauss,pars) // no extra data
{
    setParams(pars);
    setLowerBounds(lb);
    setUpperBounds(ub);
}


Gauss_Model::Gauss_Model(QVector<double> &pars): Gauss_Model(pars,emptyVector,emptyVector)
{

}


Gauss_Model::Gauss_Model(): Gauss_Model(emptyVector,emptyVector,emptyVector)
{

}

void Gauss_Model::setParams(QVector<double> &pars)
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

void Gauss_Model::setLowerBounds(QVector<double> &lb)
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

    ensureConstrains();
}


void Gauss_Model::setUpperBounds(QVector<double> &ub)
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

    ensureConstrains();
}

            /*  Gauss2D_Model class realization  */

Gauss2D_Model::Gauss2D_Model(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub, extraData_t *data):
    PSF_Model(gauss2d,pars,data)
{
    setParams(pars);
    setLowerBounds(lb);
    setUpperBounds(ub);
}


Gauss2D_Model::Gauss2D_Model(QVector<double> &pars, extraData_t *data): Gauss2D_Model(pars,emptyVector,emptyVector,data)
{

}


Gauss2D_Model::Gauss2D_Model(): Gauss2D_Model(emptyVector,emptyVector,emptyVector,nullptr)
{

}

void Gauss2D_Model::setParams(QVector<double> &pars)
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

void Gauss2D_Model::setLowerBounds(QVector<double> &lb)
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

    ensureConstrains();
}


void Gauss2D_Model::setUpperBounds(QVector<double> &ub)
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

    ensureConstrains();
}



                    /*  Moffat_Model class realization  */

Moffat_Model::Moffat_Model(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub):
    PSF_Model(moffat,pars,nullptr)
{
    setParams(pars);
    setLowerBounds(lb);
    setUpperBounds(ub);
}


Moffat_Model::Moffat_Model(QVector<double> &pars): Moffat_Model(pars,emptyVector,emptyVector)
{

}


Moffat_Model::Moffat_Model(): Moffat_Model(emptyVector,emptyVector,emptyVector)
{

}

void Moffat_Model::setParams(QVector<double> &pars)
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

void Moffat_Model::setLowerBounds(QVector<double> &lb)
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

    ensureConstrains();
}


void Moffat_Model::setUpperBounds(QVector<double> &ub)
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

    ensureConstrains();
}

            /*  Moffat2D_Model class realization  */

Moffat2D_Model::Moffat2D_Model(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub, extraData_t *data):
    PSF_Model(moffat2d,pars,data)
{
    setParams(pars);
    setLowerBounds(lb);
    setUpperBounds(ub);
}


Moffat2D_Model::Moffat2D_Model(QVector<double> &pars, extraData_t *data): Moffat2D_Model(pars,emptyVector,emptyVector,data)
{

}


Moffat2D_Model::Moffat2D_Model(): Moffat2D_Model(emptyVector,emptyVector,emptyVector,nullptr)
{

}

void Moffat2D_Model::setParams(QVector<double> &pars)
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

void Moffat2D_Model::setLowerBounds(QVector<double> &lb)
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

    ensureConstrains();
}


void Moffat2D_Model::setUpperBounds(QVector<double> &ub)
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

    ensureConstrains();
}


//            /*  PSF_Model realization   */

//PSF_Model::PSF_Model(QObject *parent) : QObject(parent)
//{

//}
