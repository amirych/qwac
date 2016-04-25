#ifndef PSF_MODEL_H
#define PSF_MODEL_H

#include "focuswidget_global.h"

#include <memory>
#include <QVector>


#define PSF_MODEL_DEFAULT_ITMAX 100



        /*     Extra data type definition   */

struct extraData_t
{
    QVector<double> x;
    QVector<double> y;

    void *extra_params;
};


            /*     PSF model function definition    */
//
//  params - PSF model function parameters (given by user)
//  func - model function values (computed by model function)
//  extra_data - extra parameters of the model function
//
typedef int (*psfModelFunc_t)(QVector<double> &params, QVector<double> &func, extraData_t *extra_data);


        /*     Base class for PSF model parameters     */

class PSF_Model
{
public:
    PSF_Model(psfModelFunc_t func);
    PSF_Model(psfModelFunc_t func, QVector<double> &pars, extraData_t *extra_data = nullptr);

    ~PSF_Model();

    virtual void setParams(QVector<double> &pars, extraData_t *extra_data = nullptr);
    virtual void setLowerBounds(QVector<double> &lb);
    virtual void setUpperBounds(QVector<double> &ub);
    void setContrains(QVector<double> &lb, QVector<double> &ub);

    QVector<double> getParams() const;
    QVector<double> getConstrains() const;
    QVector<double> getLowerBounds() const;
    QVector<double> getUpperBounds() const;

    void setArgument(QVector<double> &x); // independent variable (to compute model)
    void setArgument(QVector<double> &x, QVector<double> &y); // independent variable (to compute 2D-model)

    QVector<double> operator()();
    QVector<double> operator()(QVector<double> &x);
    QVector<double> operator()(QVector<double> &x, QVector<double> &y);

protected:
    psfModelFunc_t modelFunc;
    extraData_t* modelFuncExtraData;

    QVector<double> params, lowerBounds, upperBounds;
    QVector<double> modelFuncValue;

    void ensureConstrains();
    void objetive_function(double *pars, double *func, int n_pars, int n_func, void* data);
    virtual void compute();
};


            /*     Gauss-family PSF model parameters classes definition    */

// no extra data! order of polynom is determinated by pars length!
// order =  pars.size() - 3, order == 0 means that no background polynom is computed
class FOCUSWIDGETSHARED_EXPORT Gauss_Model: public PSF_Model
{
public:
    explicit Gauss_Model();
    Gauss_Model(QVector<double> &pars);
    Gauss_Model(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub);

    void setParams(QVector<double> &pars);
    void setLowerBounds(QVector<double> &lb);
    void setUpperBounds(QVector<double> &ub);
};



// if data = nullptr then model consists of only 2D-gaussian (without any polynominal background)
class FOCUSWIDGETSHARED_EXPORT Gauss2D_Model: public PSF_Model
{
public:
    explicit Gauss2D_Model();
    Gauss2D_Model(QVector<double> &pars, extraData_t *data = nullptr);
    Gauss2D_Model(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub, extraData_t *data = nullptr);

    void setParams(QVector<double> &pars);
    void setLowerBounds(QVector<double> &lb);
    void setUpperBounds(QVector<double> &ub);
};



        /*     Moffat-family PSF model parameters classes definition    */

// see Gauss_Model class
class FOCUSWIDGETSHARED_EXPORT Moffat_Model: public PSF_Model
{
public:
    explicit Moffat_Model();
    Moffat_Model(QVector<double> &pars);
    Moffat_Model(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub);

    void setParams(QVector<double> &pars);
    void setLowerBounds(QVector<double> &lb);
    void setUpperBounds(QVector<double> &ub);
};



class FOCUSWIDGETSHARED_EXPORT Moffat2D_Model: public PSF_Model
{
public:
    Moffat2D_Model();
    Moffat2D_Model(QVector<double> &pars, extraData_t *data = nullptr);
    Moffat2D_Model(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub, extraData_t *data = nullptr);

    void setParams(QVector<double> &pars);
    void setLowerBounds(QVector<double> &lb);
    void setUpperBounds(QVector<double> &ub);
};


#endif // PSF_MODEL_H
