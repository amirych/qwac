#ifndef PSF_MODEL_H
#define PSF_MODEL_H

#include "focuswidget_global.h"

#include <QObject>
#include <QVector>


#define PSF_MODEL_DEFAULT_ITMAX 100

// number of model parameters without possible background polynom parameters
#define PSF_MODEL_GAUSS_MIN_NPARS 3      // [Amplitude, center, fwhm]
#define PSF_MODEL_MOFFAT_MIN_NPARS 4     // [Amplitude, center, fwhm, alpha]
#define PSF_MODEL_GAUSS2D_MIN_NPARS 6    // [Amplitude, centerX, centerY, fwhmX, fwhmY, phi]
#define PSF_MODEL_MOFFAT2D_MIN_NPARS 7   // [Amplitude, centerX, centerY, fwhmX, fwhmY, alpha, phi]



        /*     Base class for PSF model parameters     */

class PSF_Model_Params
{
public:
    PSF_Model_Params();
    PSF_Model_Params(QVector<double> &pars);
    PSF_Model_Params(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub);

    virtual void setParams(QVector<double> &pars);
    virtual void setLowerBounds(QVector<double> &lb);
    virtual void setUpperBounds(QVector<double> &ub);
    void setContrains(QVector<double> &lb, QVector<double> &ub);

    QVector<double> getParams() const;
    QVector<double> getConstrains() const;
    QVector<double> getLowerBounds() const;
    QVector<double> getUpperBounds() const;


protected:
    QVector<double> params, lowerBounds, upperBounds;
};


            /*     Gauss-family PSF model parameters classes definition    */

class FOCUSWIDGETSHARED_EXPORT Gauss_Model_Params: public PSF_Model_Params
{
public:
    Gauss_Model_Params();
    Gauss_Model_Params(QVector<double> &pars);
    Gauss_Model_Params(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub);

    void setParams(QVector<double> &pars);
    void setLowerBounds(QVector<double> &lb);
    void setUpperBounds(QVector<double> &ub);
};



class FOCUSWIDGETSHARED_EXPORT Gauss2D_Model_Params: public PSF_Model_Params
{
public:
    Gauss2D_Model_Params();
    Gauss2D_Model_Params(QVector<double> &pars);
    Gauss2D_Model_Params(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub);

    void setParams(QVector<double> &pars);
    void setLowerBounds(QVector<double> &lb);
    void setUpperBounds(QVector<double> &ub);
};



        /*     Moffat-family PSF model parameters classes definition    */

class FOCUSWIDGETSHARED_EXPORT Moffat_Model_Params: public PSF_Model_Params
{
public:
    Moffat_Model_Params();
    Moffat_Model_Params(QVector<double> &pars);
    Moffat_Model_Params(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub);

    void setParams(QVector<double> &pars);
    void setLowerBounds(QVector<double> &lb);
    void setUpperBounds(QVector<double> &ub);
};



class FOCUSWIDGETSHARED_EXPORT Moffat2D_Model_Params: public PSF_Model_Params
{
public:
    Moffat2D_Model_Params();
    Moffat2D_Model_Params(QVector<double> &pars);
    Moffat2D_Model_Params(QVector<double> &pars, QVector<double> &lb, QVector<double> &ub);

    void setParams(QVector<double> &pars);
    void setLowerBounds(QVector<double> &lb);
    void setUpperBounds(QVector<double> &ub);
};



class FOCUSWIDGETSHARED_EXPORT PSF_Model : public QObject
{
    Q_OBJECT
public:
    enum PSF_MODEL_TYPE {Gauss, Gauss2D, Moffat, Moffat2D};

    explicit PSF_Model(QObject *parent = nullptr);
    PSF_Model(PSF_Model::PSF_MODEL_TYPE type, QObject *parent = nullptr);

    PSF_MODEL_TYPE getModelType() const;
    void setMaxIter(int itmax);
    int getMaxIter() const;

protected:
    PSF_MODEL_TYPE modelType;
    int maxIter; // maximum number of iterations
};

#endif // PSF_MODEL_H
