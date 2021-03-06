#include "modelfunction.h"

#include <limits>
#include <cmath>
#include <cstring>

// default constrains
#define MODELFUNCTION_BOUND_DEFAULT_MIN -std::numeric_limits<double>::infinity()
#define MODELFUNCTION_BOUND_DEFAULT_MAX std::numeric_limits<double>::infinity()


// Gaussian and Mopffat functions defaults
#define MODELFUNCTION_FWHM_DEFAULT_MIN 2.0
#define MODELFUNCTION_FWHM_DEFAULT_MAX std::numeric_limits<double>::infinity()

#define MODELFUNCTION_CENTER_DEFAULT_MIN -std::numeric_limits<double>::infinity()
#define MODELFUNCTION_CENTER_DEFAULT_MAX std::numeric_limits<double>::infinity()

// amplitude of PSF (maximal value)
#define MODELFUNCTION_AMP_DEFAULT_MIN -std::numeric_limits<double>::infinity()
#define MODELFUNCTION_AMP_DEFAULT_MAX std::numeric_limits<double>::infinity()

// Moffat's alpha-parameter
#define MODELFUNCTION_ALPHA_DEFAULT_MIN 1.0
#define MODELFUNCTION_ALPHA_DEFAULT_MAX std::numeric_limits<double>::infinity()

// model function rotation
#define MODELFUNCTION_PHI_DEFAULT_MIN 0.0
#define MODELFUNCTION_PHI_DEFAULT_MAX 0.0001
//#define MODELFUNCTION_PHI_DEFAULT_MIN -std::numeric_limits<double>::infinity()
//#define MODELFUNCTION_PHI_DEFAULT_MAX std::numeric_limits<double>::infinity()

// polynom function coefficient default constrains
#define MODELFUNCTION_POLY_DEFAULT_MIN -std::numeric_limits<double>::infinity()
#define MODELFUNCTION_POLY_DEFAULT_MAX std::numeric_limits<double>::infinity()


#define MODELFUNCTION_FWHM_DEFAULT 5.0
#define MODELFUNCTION_AMP_DEFAULT 0.0
#define MODELFUNCTION_CENTER_DEFAULT 0.0
#define MODELFUNCTION_ALPHA_DEFAULT 3.5
#define MODELFUNCTION_PHI_DEFAULT 0.0


static std::vector<double> empty_vec;

#include <QVector>
#include <QDebug>


static void generate_vector(std::vector<double> &vec, double min_val, double max_val, double step = 1.0)
{
    if ( (max_val < min_val) && step > 0.0 ) return;
    if ( (max_val > min_val) && step < 0.0 ) return;

    vec.clear();
    double val;
    for ( val = min_val; val <= max_val; val+=step ) {
        vec.push_back(val);
    }
}

static void generate_vectors(std::vector<double> &vecX, std::vector<double> &vecY, double xmin, double xmax,
                             double ymin, double ymax, double xstep = 1.0, double ystep = 1.0)
{
    if ( (xmax < xmin) && xstep > 0.0 ) return;
    if ( (xmax > xmin) && xstep < 0.0 ) return;
    if ( (ymax < ymin) && ystep > 0.0 ) return;
    if ( (ymax > ymin) && ystep < 0.0 ) return;

    vecX.clear();
    vecY.clear();

    double valX, valY;
    for ( valY = ymin; valY <= ymax; valY+=ystep ) {
        for ( valX = xmin; valX <= xmax; valX+=xstep ) {
            vecX.push_back(valX);
            vecY.push_back(valY);
        }
    }
}

int mpfit_wrapper(int m, int n, double *p, double *dy, double **dvec, void *extra) {
    AbstractModelFunction* obj = static_cast<AbstractModelFunction*>(extra);

    if ( p != obj->params.data() ) {
        for ( int i = 0; i < n; ++i ) obj->params.data()[i] = p[i];
    }

    obj->compute();

    double *ff = obj->functionValue.data();
    double *mm = obj->measurementData.data();

//    qDebug() << "MPFIT WRAPPER: par pointer = " << p << ", should be: " << obj->params.data();
//    qDebug() << "MPFIT WRAPPER: n = " << n << ", should be: " << obj->params.size();
//    for ( int i = 0; i < n; ++i ) qDebug() << "PARS: " << p[i];
    for ( size_t i = 0; i < m; ++i ) {
        dy[i] = (mm[i] - ff[i]);
//        if ( i == 100 || i == 200) qDebug() << "dy = " << dy[i];
    }

//    dvec = nullptr;

    return 0;
}


                    /* SOME OFTEN USED MODEL FUNCTIONS DEFINITIONS */


// one-dimensional polynomial function:
//   p(x) = coeffs[0] + coeffs[1]*x + coeffs[2]*x^2 + ... + coeffs[n_coeffs-1]*x^(n_coeffs-1)

static void _poly_func(double *x, double *coeffs, size_t n_x, size_t n_coeffs, double *val){
    double arg;

    for ( size_t i = 0; i < n_x; ++i ) {
        val[i] = coeffs[0];
        arg = x[i];
        for ( size_t ord = 1; ord < n_coeffs; ++ord, arg *= arg) {
            val[i] += coeffs[ord]*arg;
        }
    }
}

// two-dimensional polynomial function:
// p(x,y) = coefss[0] + coeffs[1]*y + coeffs[2]*y^2 + ... + coeffs[yd]*y^yd +
//        + coeffs[yd+1]*x + coeffs[yd+2]*x*y + coeffs[yd+3]*x*y^2 + ... + coeffs[2*(yd+1)-1]*x*y^yd +
//        + coeffs[2*(yd+1)]*x^2 + coeffs[2*(yd+1)+1]*x^2*y + coeffs[2*(yd+1)+2]*x^2*y^2 + ... + coeffs[3*(yd+1)-1]*x^2*y^yd +
//        ...
//        + coeffs[xd*(yd+1)]*x^xd + coeffs[xd*(yd+1)+1]*x^xd*y + coeffs[xd*(yd+1)+2]*x^xd*y^2 + ... + coeffs[(xd+1)*(y+1)-1]*x^xd*y^yd
//
//  , where xd and yd are polynom degree along X and Y

static void _poly2d_func(double *x, double *y, double *coeffs, size_t n, size_t x_degree, size_t y_degree, double *val){
    double argx,argy;
    size_t idx;

    for ( size_t i = 0; i < n; ++i ) {
//        val[i] = coeffs[0];
        val[i] = 0.0;
        idx = 0;
        argx = 1.0;
        for ( size_t ord_x = 0; ord_x <= x_degree; ++ord_x, argx *= x[i] ) {
            argy = 1.0;
            for ( size_t ord_y = 0; ord_y <= y_degree; ++ord_y, argy *= y[i] ) {
                val[i] += coeffs[idx]*argx*argy;
//                if ( i == 100) qDebug() << "poly2d: " << val[i];
//                if ( !std::isfinite(val[i]) ) qDebug() << "poly2d: bad value at " << i << ", ord_x = " << ord_x << ", ord_y = " << ord_y << ", coeffs = " << coeffs[idx];
                ++idx;
            }

        }
    }
//    qDebug() << "poly2d coeff[0] = " << coeffs[0];
}


static int poly_func(std::vector<double> &x, std::vector<double> &coeffs, std::vector<double> &func, void *extra_pars)
{
    _poly_func(x.data(), coeffs.data(), x.size(), coeffs.size(), func.data());

    return 0;
}


static int poly2d_func(std::vector<double> &x, std::vector<double> &y, std::vector<double> &coeffs, std::vector<double> &func, void *extra_pars)
{
    if ( extra_pars == nullptr ) return 1;
    size_t* deg = (size_t*)extra_pars; // must be at least 2-element size_t array
    if ( coeffs.size() < ((deg[0]+1)*(deg[1]+1))) return 2; // not enough coefficients

    _poly2d_func(x.data(), y.data(), coeffs.data(), x.size(), deg[0], deg[1], func.data());

    return 0;
}


//
//  One-dimensional Gaussian:
//    g(x) = a0*exp(-z^2/2.0), where z = (x-a1)/sigma and sigma = a2/2.0/sqrt(2.0*ln(2.0))
//    a0 - amplitude
//    a1 - center
//    a2 - FWHM
//
//  It is possible to compute Gaussian with optional polynomial background.
//  If number of pars is greater than MODELFUNCTION_GAUSS_MIN_NPARS then
//  extra pars are interpreted as coefficients of polynom (see poly_func function above)

static int gaussian_func(std::vector<double> &x, std::vector<double> &pars, std::vector<double> &func, void *extra_pars)
{
    double sigma = pars[2]/2.35482004503; // FWHM to gaussian sigma

    // compute background polynom
    if ( pars.size() > MODELFUNCTION_GAUSS_MIN_NPARS ) {
        size_t n_poly = pars.size() - MODELFUNCTION_GAUSS_MIN_NPARS;
        _poly_func(x.data(),pars.data()+MODELFUNCTION_GAUSS_MIN_NPARS,x.size(),n_poly,func.data());
    } else {
        memset(func.data(),0,func.size()*sizeof(double));
    }

    // compute Gaussian
    for ( size_t i = 0; i < x.size(); ++i ) {
        double z2 = (x[i] - pars[1])/sigma;
        z2 *= z2;
        func[i] += pars[0]*exp(-z2/2.0);
    }

    return 0;
}


//
//  Moffat one-dimensional function:
//    M(x) = a0*(1+z^2)^(-a4), where z = (x-a1)/s and s = a2/2/sqrt(2^(1/a4)-1)
//    a0 - amplitude
//    a1 - center
//    a2 - FWHM
//    a3 - moffat's parameter
//
//  It is possible to compute Moffat function with optional polynomial background.
//  If number of pars is greater than MODELFUNCTION_MOFFAT_MIN_NPARS then
//  extra pars are interpreted as coefficients of polynom (see poly_func function above)

static int moffat_func(std::vector<double> &x, std::vector<double> &pars, std::vector<double> &func, void *extra_pars)
{
    if ( pars.size() > MODELFUNCTION_MOFFAT_MIN_NPARS ) {
        size_t n_poly = pars.size() - MODELFUNCTION_MOFFAT_MIN_NPARS;
        _poly_func(x.data(),pars.data()+MODELFUNCTION_MOFFAT_MIN_NPARS,x.size(),n_poly,func.data());
    } else {
        memset(func.data(),0,func.size()*sizeof(double));
    }

    // compute Moffat function
    for ( size_t i = 0; i < x.size(); ++i ) {
        double term = pars[2]/2.0/sqrt(pow(2.0,1.0/pars[3])-1.0);
        double z2 = (x[i] - pars[1])/term;
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

     pars = [A, xc, yc, fwhmX, fwhmY, phi, ...], where ... denotes possible additive 2D-polynom coefficients (see poly2d_func)
            phi is in radians

*/

static int gaussian2d_func(std::vector<double> &x, std::vector<double> &y, std::vector<double> &pars, std::vector<double> &func, void *extra_pars)
{
    if ( pars.size() > MODELFUNCTION_GAUSS2D_MIN_NPARS ) { // compute bivariate background polynom
        if ( extra_pars == nullptr ) return 1;
        size_t* deg = (size_t*)extra_pars; // must be at least 2-element size_t array
        size_t n_poly = pars.size() - MODELFUNCTION_GAUSS2D_MIN_NPARS;
        if ( n_poly < ((deg[0]+1)*(deg[1]+1))) return 2; // not enough coefficients
        _poly2d_func(x.data(),y.data(),pars.data()+ MODELFUNCTION_GAUSS2D_MIN_NPARS,x.size(),deg[0],deg[1],func.data());
    } else {
        memset(func.data(),0,func.size()*sizeof(double));
    }

    double k = 4.0*log(2.0);

    for ( size_t i = 0; i < x.size(); ++i ) {
        // first, rotate cordinate system
        double t = x[i] - pars[1];
        double yy = y[i] - pars[2];

        double xx = t*cos(pars[5]) - yy*sin(pars[5]);
        yy = yy*cos(pars[5]) + t*sin(pars[5]);

        double r1 = xx/pars[3];
        double r2 = yy/pars[4];
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

static int moffat2d_func(std::vector<double> &x, std::vector<double> &y, std::vector<double> &pars, std::vector<double> &func, void *extra_pars)
{
    if ( pars.size() > MODELFUNCTION_MOFFAT2D_MIN_NPARS ) { // compute bivariate background polynom
        if ( extra_pars == nullptr ) return 1;
        size_t* deg = (size_t*)extra_pars; // must be at least 2-element size_t array
        size_t n_poly = pars.size() - MODELFUNCTION_MOFFAT2D_MIN_NPARS;
//        qDebug() << "n_poly = " << n_poly << ", n_pars = " << pars.size() << ", must be: " << (deg[0]+1)*(deg[1]+1);
        if ( n_poly < ((deg[0]+1)*(deg[1]+1))) return 2; // not enough coefficients
        _poly2d_func(x.data(),y.data(),pars.data()+ MODELFUNCTION_MOFFAT2D_MIN_NPARS,x.size(),deg[0],deg[1],func.data());
    } else {
        memset(func.data(),0,func.size()*sizeof(double));
    }

//    qDebug() << "poly val: " << func[100];

    double k = 2.0*(pow(2.0,1.0/pars[5])-1.0);

    for ( size_t i = 0; i < x.size(); ++i ) {
        // first, rotate cordinate system
        double t = x[i] - pars[1];
        double yy = y[i] - pars[2];

        double xx = t*cos(pars[6]) - yy*sin(pars[6]);
        yy = yy*cos(pars[6]) + t*sin(pars[6]);

        double r1 = xx/pars[3];
        double r2 = yy/pars[4];

//        double r1 = (x[i]-pars[1])/pars[3];
//        double r2 = (y[i]-pars[2])/pars[4];
        double z2 = r1*r1 + r2*r2;

        func[i] += pars[0]*pow(1.0+k*z2,-pars[5]);
//        if ( !std::isfinite(func[i]) ) qDebug() << "moffat2d: bad value at " << i;
    }

    return 0;
}


                /* BASE ABSTRACT CLASS REALIZATION */

AbstractModelFunction::AbstractModelFunction(std::vector<double> &pars, void *extra_pars):
    params(pars), extraParams(extra_pars),
    lowerBounds(std::vector<double>()), upperBounds(std::vector<double>()),
    fixedParamNumber(std::vector<size_t>()),
    functionValue(std::vector<double>()), measurementData(std::vector<double>()),
    modelFunctionName(""), maxIter(MODELFUNCTION_DEFAULT_MAXITER),
    fitStatus(0), compStatus(0)
{
    checkParams();

    memset(&fitInfo,0,sizeof(mp_result));
    memset(&fitConfig,0,sizeof(mp_config));

    fitConfig.maxiter = maxIter;
}


AbstractModelFunction::~AbstractModelFunction()
{

}

void AbstractModelFunction::setParams(std::vector<double> &pars, void *extra_pars)
{
    params = pars;
    extraParams = extra_pars;

    checkParams();
}


void AbstractModelFunction::setConstrains(std::vector<double> &lb, std::vector<double> &ub)
{
    lowerBounds = lb;
    upperBounds = ub;

    checkConstrains();
}


void AbstractModelFunction::setFixedParams(std::vector<size_t> &param_num)
{
    fixedParamNumber = param_num;
}


std::vector<double> AbstractModelFunction::getParams()
{
    return params;
}


void AbstractModelFunction::getConstrains(std::vector<double> *lb, std::vector<double> *ub)
{
    checkConstrains();
    if ( lb ) *lb = lowerBounds;
    if ( ub ) *ub = upperBounds;
}


std::vector<size_t> AbstractModelFunction::getFixedParams() const
{
    return fixedParamNumber;
}


void AbstractModelFunction::setModelFunctionName(std::string &name)
{
    modelFunctionName = name;
}


void AbstractModelFunction::setModelFunctionName(std::string name)
{
    modelFunctionName = name;
}


void AbstractModelFunction::fitData(std::vector<double> &data)
{
    measurementData = data;

    fitting();
}

void AbstractModelFunction::fitData(double *data, size_t data_len)
{
    measurementData = std::vector<double>(data, data+data_len);

    fitting();
}


void AbstractModelFunction::fitting()
{
    checkConstrains();

    mp_par pars[params.size()];
    memset(&pars,0,sizeof(pars));

    for ( size_t i = 0; i < params.size(); ++i ) {
        pars[i].fixed = 0;
        pars[i].limited[0] = 1;
        pars[i].limited[1] = 1;
        pars[i].limits[0] = lowerBounds[i];
        pars[i].limits[1] = upperBounds[i];
    }

    for ( size_t i = 0; i < fixedParamNumber.size(); ++i ) {
        size_t idx = fixedParamNumber[i];
        if ( idx < params.size() ) pars[idx].fixed = 1;
    }

    int fitStatus = mpfit(mpfit_wrapper,measurementData.size(),params.size(),params.data(),pars,&fitConfig,(void*)this,&fitInfo);
//    int fitStatus = mpfit(mpfit_wrapper,measurementData.size(),params.size(),params.data(),pars,0,(void*)this,&fitInfo);

//    qDebug() << "MPFIT RES: status = " << fitStatus << ", bestnorm: " << fitInfo.bestnorm << ", niter = " << fitInfo.niter;
//    qDebug() << "MPFIT PARS: " << QVector<double>::fromStdVector(params);

}


std::vector<double> AbstractModelFunction::operator ()()
{
    compute();

    return functionValue;
}


int AbstractModelFunction::getFitStatus() const
{
    return fitStatus;
}


void AbstractModelFunction::getFitInfo(mp_result *fit_info)
{
    if ( !fit_info ) return;
    memcpy(fit_info,&fitInfo,sizeof(mp_result));
}


void AbstractModelFunction::setMaxIterations(int max_iter)
{
    maxIter = max_iter > 0 ? max_iter : MODELFUNCTION_DEFAULT_MAXITER;
    fitConfig.maxiter = maxIter;
}


int AbstractModelFunction::getCompStatus() const
{
    return compStatus;
}


void AbstractModelFunction::checkParams()
{
}


void AbstractModelFunction::checkConstrains()
{
    size_t n = lowerBounds.size();
    if ( n < params.size() ) { // set defaults
        lowerBounds.resize(params.size());
        for ( size_t i = n; i < params.size(); ++i ) lowerBounds[i] = MODELFUNCTION_BOUND_DEFAULT_MIN;
    }

    n = upperBounds.size();
    if ( n < params.size() ) { // set defaults
        upperBounds.resize(params.size());
        for ( size_t i = n; i < params.size(); ++i ) upperBounds[i] = MODELFUNCTION_BOUND_DEFAULT_MAX;
    }

}

                /* ONE-DIMENSIONAL MODEL FUNCTION CLASS REALIZATION */


ModelFunction::ModelFunction(user_func_t user_func, std::vector<double> &pars, void *extra_pars): AbstractModelFunction(pars,extra_pars),
    userFunction(user_func)
{

}


ModelFunction::ModelFunction(user_func_t user_func): ModelFunction(user_func,empty_vec,nullptr)
{

}


ModelFunction::~ModelFunction()
{

}


void ModelFunction::setArgument(std::vector<double> &arg)
{
    argumentX = arg;
    if ( functionValue.size() < argumentX.size() ) functionValue.resize(argumentX.size());
}

void ModelFunction::setArgument(double xmin, double xmax, double xstep)
{
    generate_vector(argumentX,xmin,xmax,xstep);
    if ( functionValue.size() < argumentX.size() ) functionValue.resize(argumentX.size());
}


void ModelFunction::compute()
{
    if ( !argumentX.size() ) {
        compStatus = -1;
        return;
    }

    compStatus = userFunction(argumentX,params,functionValue,extraParams);
}


std::vector<double> ModelFunction::operator ()(std::vector<double> &x)
{
    setArgument(x);

    return AbstractModelFunction::operator()();
}

                /* TWO-DIMENSIONAL MODEL FUNCTION CLASS REALIZATION */


ModelFunction2D::ModelFunction2D(user_func2D_t user_func, std::vector<double> &pars, void *extra_pars): AbstractModelFunction(pars,extra_pars),
    userFunction(user_func)
{
}


ModelFunction2D::ModelFunction2D(user_func2D_t user_func): ModelFunction2D(user_func,empty_vec,nullptr)
{

}


ModelFunction2D::~ModelFunction2D()
{

}


void ModelFunction2D::setArgument(std::vector<double> &argX, std::vector<double> &argY)
{
    compStatus = 0;
    if ( argX.size() == argY.size() ) {
        argumentX = argX;
        argumentY = argY;
    } else {
        compStatus = -1;
        return;
    }
    if ( functionValue.size() < argumentX.size() ) functionValue.resize(argumentX.size());
}


void ModelFunction2D::setArgument(double xmin, double xmax, double ymin, double ymax, double xstep, double ystep)
{
    generate_vectors(argumentX,argumentY,xmin,xmax,ymin,ymax,xstep,ystep);

    if ( functionValue.size() < argumentX.size() ) functionValue.resize(argumentX.size());
}


void ModelFunction2D::compute()
{
    if ( !argumentX.size() || !argumentY.size() ) {
        compStatus = -1;
        return;
    }

    compStatus = userFunction(argumentX,argumentY,params,functionValue,extraParams);
}


std::vector<double> ModelFunction2D::operator ()(std::vector<double> &x, std::vector<double> &y)
{
    setArgument(x,y);

    return AbstractModelFunction::operator()();
}


                            /**********************************************/
                            /* GAUSSIAN AND MOFFAT MODEL FUNCTION CLASSES */
                            /**********************************************/


                                    /* ONE-DIMENSIONAL FUNCTIONS */

GaussModelFunction::GaussModelFunction(std::vector<double> &pars, void *extra_pars): ModelFunction(gaussian_func, pars, extra_pars)
{
    checkParams();
}


GaussModelFunction::GaussModelFunction(): GaussModelFunction(empty_vec, nullptr)
{

}


void GaussModelFunction::checkParams()
{
    if ( params.size() < MODELFUNCTION_GAUSS_MIN_NPARS ) {
        int N = params.size();
        params.resize(MODELFUNCTION_GAUSS_MIN_NPARS);
        int i = MODELFUNCTION_GAUSS_MIN_NPARS;
        if ( N < i ) params[--i] = MODELFUNCTION_FWHM_DEFAULT;
        if ( N < i ) params[--i] = MODELFUNCTION_CENTER_DEFAULT;
        if ( N == 0 ) params[0] = MODELFUNCTION_AMP_DEFAULT;
    }

    if ( extraParams != nullptr ) { // resize params vector if needed
        size_t *n_poly = (size_t*)extraParams;
        size_t len = MODELFUNCTION_GAUSS_MIN_NPARS + (*n_poly+1);
        if ( params.size() < len ) params.resize(len);
    }
}


void GaussModelFunction::checkConstrains()
{
    if ( lowerBounds.size() < MODELFUNCTION_GAUSS_MIN_NPARS ) { // set default constrains
        int N = lowerBounds.size();
        lowerBounds.resize(MODELFUNCTION_GAUSS_MIN_NPARS);
        int i = MODELFUNCTION_GAUSS_MIN_NPARS;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_FWHM_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_CENTER_DEFAULT_MIN;
        if ( N == 0 ) lowerBounds[0] = MODELFUNCTION_AMP_DEFAULT_MIN;
    }

    if ( lowerBounds.size() < params.size() ) { // set default polynomial constrains
        lowerBounds.resize(params.size());
        for ( size_t i = MODELFUNCTION_GAUSS_MIN_NPARS; i < params.size(); ++i ) lowerBounds[i] = MODELFUNCTION_POLY_DEFAULT_MIN;
    }

    if ( upperBounds.size() < MODELFUNCTION_GAUSS_MIN_NPARS ) { // set default constrains
        int N = upperBounds.size();
        upperBounds.resize(MODELFUNCTION_GAUSS_MIN_NPARS);
        int i = MODELFUNCTION_GAUSS_MIN_NPARS;
        if ( N < i ) upperBounds[--i] = MODELFUNCTION_FWHM_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = MODELFUNCTION_CENTER_DEFAULT_MAX;
        if ( N == 0 ) upperBounds[0] = MODELFUNCTION_AMP_DEFAULT_MAX;
    }

    if ( upperBounds.size() < params.size() ) { // set default polynomial constrains
        upperBounds.resize(params.size());
        for ( size_t i = MODELFUNCTION_GAUSS_MIN_NPARS; i < params.size(); ++i ) upperBounds[i] = MODELFUNCTION_POLY_DEFAULT_MAX;
    }

}


MoffatModelFunction::MoffatModelFunction(std::vector<double> &pars, void *extra_pars): ModelFunction(moffat_func, pars, extra_pars)
{
    checkParams();
}


MoffatModelFunction::MoffatModelFunction(): MoffatModelFunction(empty_vec, nullptr)
{

}


void MoffatModelFunction::checkParams()
{
    if ( params.size() < MODELFUNCTION_MOFFAT_MIN_NPARS ) {
        int N = params.size();
        params.resize(MODELFUNCTION_MOFFAT_MIN_NPARS);
        int i = MODELFUNCTION_MOFFAT_MIN_NPARS;
        if ( N < i ) params[--i] = MODELFUNCTION_ALPHA_DEFAULT;
        if ( N < i ) params[--i] = MODELFUNCTION_FWHM_DEFAULT;
        if ( N < i ) params[--i] = MODELFUNCTION_CENTER_DEFAULT;
        if ( N == 0 ) params[0] = MODELFUNCTION_AMP_DEFAULT;
    }

    if ( extraParams != nullptr ) { // resize params vector if needed
        size_t *n_poly = (size_t*)extraParams;
        size_t len = MODELFUNCTION_MOFFAT_MIN_NPARS + (*n_poly+1);
        if ( params.size() < len ) params.resize(len);
    }
}


void MoffatModelFunction::checkConstrains()
{
    if ( lowerBounds.size() < MODELFUNCTION_MOFFAT_MIN_NPARS ) { // set default constrains
        int N = lowerBounds.size();
        lowerBounds.resize(MODELFUNCTION_MOFFAT_MIN_NPARS);
        int i = MODELFUNCTION_MOFFAT_MIN_NPARS;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_ALPHA_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_FWHM_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_CENTER_DEFAULT_MIN;
        if ( N == 0 ) lowerBounds[0] = MODELFUNCTION_AMP_DEFAULT_MIN;
    }

    if ( lowerBounds.size() < params.size() ) { // set default polynomial constrains
        lowerBounds.resize(params.size());
        for ( size_t i = MODELFUNCTION_MOFFAT_MIN_NPARS; i < params.size(); ++i ) lowerBounds[i] = MODELFUNCTION_POLY_DEFAULT_MIN;
    }

    if ( upperBounds.size() < MODELFUNCTION_MOFFAT_MIN_NPARS ) { // set default constrains
        int N = upperBounds.size();
        upperBounds.resize(MODELFUNCTION_MOFFAT_MIN_NPARS);
        int i = MODELFUNCTION_MOFFAT_MIN_NPARS;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_ALPHA_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = MODELFUNCTION_FWHM_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = MODELFUNCTION_CENTER_DEFAULT_MAX;
        if ( N == 0 ) upperBounds[0] = MODELFUNCTION_AMP_DEFAULT_MAX;
    }

    if ( upperBounds.size() < params.size() ) { // set default polynomial constrains
        upperBounds.resize(params.size());
        for ( size_t i = MODELFUNCTION_MOFFAT_MIN_NPARS; i < params.size(); ++i ) upperBounds[i] = MODELFUNCTION_POLY_DEFAULT_MAX;
    }

}


                                /* TWO-DIMENSIONAL FUNCTIONS */

Gauss2DModelFunction::Gauss2DModelFunction(std::vector<double> &pars, void *extra_pars): ModelFunction2D(gaussian2d_func, pars, extra_pars)
{
    checkParams();
    modelFunctionName = "Guass2D";
}


Gauss2DModelFunction::Gauss2DModelFunction(): Gauss2DModelFunction(empty_vec, nullptr)
{

}

void Gauss2DModelFunction::checkParams()
{
    if ( params.size() < MODELFUNCTION_GAUSS2D_MIN_NPARS ) {
        int N = params.size();
        params.resize(MODELFUNCTION_GAUSS2D_MIN_NPARS);
        int i = MODELFUNCTION_GAUSS2D_MIN_NPARS;
        if ( N < i ) params[--i] = MODELFUNCTION_PHI_DEFAULT;
        if ( N < i ) params[--i] = MODELFUNCTION_FWHM_DEFAULT;
        if ( N < i ) params[--i] = MODELFUNCTION_FWHM_DEFAULT;
        if ( N < i ) params[--i] = MODELFUNCTION_CENTER_DEFAULT;
        if ( N < i ) params[--i] = MODELFUNCTION_CENTER_DEFAULT;
        if ( N == 0 ) params[0] = MODELFUNCTION_AMP_DEFAULT;
    }

    if ( extraParams != nullptr ) { // resize params vector if needed
        size_t *deg = (size_t*)extraParams;
        size_t n_poly = (deg[0]+1)*(deg[1]+1);
        size_t len = MODELFUNCTION_GAUSS2D_MIN_NPARS + n_poly;
        if ( params.size() < len ) params.resize(len);
    }
}


void Gauss2DModelFunction::checkConstrains()
{
    if ( lowerBounds.size() < MODELFUNCTION_GAUSS2D_MIN_NPARS ) { // set default constrains
        int N = lowerBounds.size();
        lowerBounds.resize(MODELFUNCTION_GAUSS2D_MIN_NPARS);
        int i = MODELFUNCTION_GAUSS2D_MIN_NPARS;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_PHI_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_FWHM_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_FWHM_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_CENTER_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_CENTER_DEFAULT_MIN;
        if ( N == 0 ) lowerBounds[0] = MODELFUNCTION_AMP_DEFAULT_MIN;
    }

    if ( lowerBounds.size() < params.size() ) { // set default polynomial constrains
        lowerBounds.resize(params.size());
        for ( size_t i = MODELFUNCTION_GAUSS2D_MIN_NPARS; i < params.size(); ++i ) lowerBounds[i] = MODELFUNCTION_POLY_DEFAULT_MIN;
    }

    if ( upperBounds.size() < MODELFUNCTION_GAUSS2D_MIN_NPARS ) { // set default constrains
        int N = upperBounds.size();
        upperBounds.resize(MODELFUNCTION_GAUSS2D_MIN_NPARS);
        int i = MODELFUNCTION_GAUSS2D_MIN_NPARS;
        if ( N < i ) upperBounds[--i] = MODELFUNCTION_PHI_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = MODELFUNCTION_FWHM_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = MODELFUNCTION_FWHM_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = MODELFUNCTION_CENTER_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = MODELFUNCTION_CENTER_DEFAULT_MAX;
        if ( N == 0 ) upperBounds[0] = MODELFUNCTION_AMP_DEFAULT_MAX;
    }

    if ( upperBounds.size() < params.size() ) { // set default polynomial constrains
        upperBounds.resize(params.size());
        for ( size_t i = MODELFUNCTION_GAUSS2D_MIN_NPARS; i < params.size(); ++i ) upperBounds[i] = MODELFUNCTION_POLY_DEFAULT_MAX;
    }

}




Moffat2DModelFunction::Moffat2DModelFunction(std::vector<double> &pars, void *extra_pars): ModelFunction2D(moffat2d_func, pars, extra_pars)
{
    checkParams();
    modelFunctionName = "Moffat2D";
}


Moffat2DModelFunction::Moffat2DModelFunction(): Moffat2DModelFunction(empty_vec, nullptr)
{

}


void Moffat2DModelFunction::checkParams()
{
    if ( params.size() < MODELFUNCTION_MOFFAT2D_MIN_NPARS ) {
        int N = params.size();
        params.resize(MODELFUNCTION_MOFFAT2D_MIN_NPARS);
        int i = MODELFUNCTION_MOFFAT2D_MIN_NPARS;
        if ( N < i ) params[--i] = MODELFUNCTION_PHI_DEFAULT;
        if ( N < i ) params[--i] = MODELFUNCTION_ALPHA_DEFAULT;
        if ( N < i ) params[--i] = MODELFUNCTION_FWHM_DEFAULT;
        if ( N < i ) params[--i] = MODELFUNCTION_FWHM_DEFAULT;
        if ( N < i ) params[--i] = MODELFUNCTION_CENTER_DEFAULT;
        if ( N < i ) params[--i] = MODELFUNCTION_CENTER_DEFAULT;
        if ( N == 0 ) params[0] = MODELFUNCTION_AMP_DEFAULT;
    }

    if ( extraParams != nullptr ) { // resize params vector if needed
        size_t *deg = (size_t*)extraParams;
        size_t n_poly = (deg[0]+1)*(deg[1]+1);
        size_t len = MODELFUNCTION_MOFFAT2D_MIN_NPARS + n_poly;
        if ( params.size() < len ) params.resize(len);
    }

}


void Moffat2DModelFunction::checkConstrains()
{
    if ( lowerBounds.size() < MODELFUNCTION_MOFFAT2D_MIN_NPARS ) { // set default constrains
        int N = lowerBounds.size();
        lowerBounds.resize(MODELFUNCTION_MOFFAT2D_MIN_NPARS);
        int i = MODELFUNCTION_MOFFAT2D_MIN_NPARS;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_PHI_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_ALPHA_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_FWHM_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_FWHM_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_CENTER_DEFAULT_MIN;
        if ( N < i ) lowerBounds[--i] = MODELFUNCTION_CENTER_DEFAULT_MIN;
        if ( N == 0 ) lowerBounds[0] = MODELFUNCTION_AMP_DEFAULT_MIN;
    }

    if ( lowerBounds.size() < params.size() ) { // set default polynomial constrains
        lowerBounds.resize(params.size());
        for ( size_t i = MODELFUNCTION_MOFFAT2D_MIN_NPARS; i < params.size(); ++i ) lowerBounds[i] = MODELFUNCTION_POLY_DEFAULT_MIN;
    }

    if ( upperBounds.size() < MODELFUNCTION_MOFFAT2D_MIN_NPARS ) { // set default constrains
        int N = upperBounds.size();
        upperBounds.resize(MODELFUNCTION_MOFFAT2D_MIN_NPARS);
        int i = MODELFUNCTION_MOFFAT2D_MIN_NPARS;
        if ( N < i ) upperBounds[--i] = MODELFUNCTION_PHI_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = MODELFUNCTION_ALPHA_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = MODELFUNCTION_FWHM_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = MODELFUNCTION_FWHM_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = MODELFUNCTION_CENTER_DEFAULT_MAX;
        if ( N < i ) upperBounds[--i] = MODELFUNCTION_CENTER_DEFAULT_MAX;
        if ( N == 0 ) upperBounds[0] = MODELFUNCTION_AMP_DEFAULT_MAX;
    }

    if ( upperBounds.size() < params.size() ) { // set default polynomial constrains
        upperBounds.resize(params.size());
        for ( size_t i = MODELFUNCTION_MOFFAT2D_MIN_NPARS; i < params.size(); ++i ) upperBounds[i] = MODELFUNCTION_POLY_DEFAULT_MAX;
    }
}


                    /* Polynomial functions classes */

Polynomial::Polynomial(std::vector<double> &coeffs): ModelFunction(poly_func,coeffs)
{
    modelFunctionName = "Polynomial";
}


Polynomial::Polynomial(): Polynomial(empty_vec)
{

}


void Polynomial::compute()
{
    if ( params.size() == 0 ) {
        functionValue.assign(functionValue.size(),0.0);
        return;
    }

    ModelFunction::compute();
}


void Polynomial::checkConstrains()
{
    if ( lowerBounds.size() < params.size() ) {
        size_t N = lowerBounds.size();
        lowerBounds.resize(params.size());
        for (size_t i = N; i < params.size(); ++i ) lowerBounds[i] = MODELFUNCTION_POLY_DEFAULT_MIN;
    }

    if ( upperBounds.size() < params.size() ) {
        size_t N = upperBounds.size();
        upperBounds.resize(params.size());
        for (size_t i = N; i < params.size(); ++i ) upperBounds[i] = MODELFUNCTION_POLY_DEFAULT_MAX;
    }
}



Polynomial2D::Polynomial2D(std::vector<double> &coeffs, std::vector<size_t> &degree): ModelFunction2D(poly2d_func, coeffs)
{
    if ( degree.size() >= 2 ) {
        polyDegree[0] = degree[0];
        polyDegree[1] = degree[1];
    } else if ( degree.size() == 1 ) {
        polyDegree[0] = degree[0];
        polyDegree[1] = degree[0];
    } else {
        polyDegree[0] = 0;
        polyDegree[1] = 0;
    }

    extraParams = polyDegree;
    modelFunctionName = "Polynomial2D";
}


static std::vector<size_t> empty_vec_int;
Polynomial2D::Polynomial2D(): Polynomial2D(empty_vec,empty_vec_int)
{

}


void Polynomial2D::compute()
{
    if ( params.size() == 0 ) {
        functionValue.assign(functionValue.size(),0.0);
        return;
    }

    ModelFunction2D::compute();
}
