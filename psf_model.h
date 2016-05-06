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
//   The type definition follows "levmar" package definition
//   for objective function. According this, the last parameter
//   type is void*, but PSF_Model class realization assumes that
//   user-written function uses of extraData_t structure pointer.
//
typedef void (*psfModelFunc_t)(double* params, double* func, int n_params, int n_func, void *extra_data);


        /*     Base class for PSF model parameters     */

class PSF_Model
{
public:
    PSF_Model(psfModelFunc_t func);
    PSF_Model(psfModelFunc_t func, QVector<double> &pars, extraData_t *extra_data = nullptr);

//    ~PSF_Model();

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

    QVector<double> operator()(); // return computed model values
    QVector<double> operator()(QVector<double> &x);
    QVector<double> operator()(QVector<double> &x, QVector<double> &y);

    void fitData(QVector<double> &data); // fit model to data
    void fitData(QVector<double> &x, QVector<double> &data);
    void fitData(QVector<double> &x, QVector<double> &y, QVector<double> &data);

    void setMaxIter(int max_iter);
    int getMaxIter() const;

    void setPsfModelName(QString &name);
    QString getPsfModelName() const;

//    void objective_function(double *pars, double *func, int n_pars, int n_func, void* data);

protected:
    QString modelName;
    psfModelFunc_t modelFunc;
    extraData_t* modelFuncExtraData;

    QVector<double> params, lowerBounds, upperBounds;
    QVector<double> modelFuncValue;

    int maxIter;

    void ensureConstrains();
    virtual void compute();
};

//std::function<void(double*, double*, int, int, void*)>
//  modelFuncCall(PSF_Model& obj)
//{
//  return [&](double* a, double* b, int c, int d, void* e){
//    return obj.objective_function(a, b, c, d, e); };
//}

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
