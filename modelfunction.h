#ifndef MODELFUNCTION_H
#define MODELFUNCTION_H

#include <mpfit.h>
#include <vector>
#include <string>


                            /* DEFAULT DEFINITIONS */

#define MODELFUNCTION_DEFAULT_MAXITER 300

// number of model parameters for predefined functions without possible background polynom parameters
#define MODELFUNCTION_GAUSS_MIN_NPARS 3      // [Amplitude, center, fwhm]
#define MODELFUNCTION_MOFFAT_MIN_NPARS 4     // [Amplitude, center, fwhm, alpha]
#define MODELFUNCTION_GAUSS2D_MIN_NPARS 6    // [Amplitude, centerX, centerY, fwhmX, fwhmY, phi]
#define MODELFUNCTION_MOFFAT2D_MIN_NPARS 7   // [Amplitude, centerX, centerY, fwhmX, fwhmY, alpha, phi]

                        /* USER FUNCTION TYPES DECLARATION */

/*
 NOTES:

 The user function must return 0 if computations is OK or >0 if error occured
 returned values <0 are reserved by model function classes:
   -1 : invalid arguments (length of argument vectors is 0 or is not equal (for 2D functions) )


 It is guaranteed that func-vector passed to user function has the same length as for arguments vectors!
 (i.e. it should not verified inside user-written function)

*/


// one-dimension function
typedef int (*user_func_t)(std::vector<double> &x, std::vector<double> &pars, std::vector<double> &func, void *extra_pars);

// two-dimension function
typedef int (*user_func2D_t)(std::vector<double> &x, std::vector<double> &y, std::vector<double> &pars, std::vector<double> &func, void *extra_pars);


                            /* BASE ABSTRACT CLASS */

//void levmar_obj_func_wrapper(double *pars, double*func, int n_pars, int n_func, void* extra);

class AbstractModelFunction
{
public:
    friend void levmar_obj_func_wrapper(double *pars, double*func, int n_pars, int n_func, void* extra);
    friend int mpfit_wrapper(int m, int n, double *p, double *dy, double **dvec, void *extra);

    AbstractModelFunction(std::vector<double> &pars, void *extra_pars = nullptr);
    virtual ~AbstractModelFunction();

    virtual void setParams(std::vector<double> &pars, void *extra_pars = nullptr);
    virtual std::vector<double> getParams();
    virtual void setConstrains(std::vector<double> &lb, std::vector<double> &ub);
    virtual void getConstrains(std::vector<double> *lb, std::vector<double> *ub);
    virtual void setFixedParams(std::vector<size_t> &param_num);
    std::vector<size_t> getFixedParams() const;

    void setModelFunctionName(std::string &name);
    void setModelFunctionName(std::string name);

    void fitData(std::vector<double> &data);
    void fitData(double *data, size_t data_len);

    std::vector<double> operator()();

    int getFitStatus() const;
    void getFitInfo(mp_result *fit_info);
    void setMaxIterations(int max_iter);

    int getCompStatus() const;

protected:
    std::vector<double> params;
    void *extraParams;
    std::vector<double> lowerBounds, upperBounds;
    std::vector<size_t> fixedParamNumber;
    std::vector<double> functionValue;
    std::vector<double> measurementData;

    std::string modelFunctionName;

    int maxIter;

    int fitStatus;
    mp_result fitInfo;
    mp_config fitConfig;

    int compStatus; // status of user function computations

    virtual void compute() = 0;

    virtual void checkParams();
    virtual void checkConstrains();

    virtual void fitting();
};


                    /* ONE-DIMENSIONAL MODEL FUNCTION BASE CLASS */

class ModelFunction: public AbstractModelFunction
{
public:
    ModelFunction(user_func_t user_func);
    ModelFunction(user_func_t user_func, std::vector<double> &pars, void *extra_pars = nullptr);

    virtual ~ModelFunction();

    void setArgument(std::vector<double> &arg);
    void setArgument(double xmin, double xmax, double xstep = 1.0);

    std::vector<double> operator()(std::vector<double> &x);

protected:
    user_func_t userFunction;
    std::vector<double> argumentX;

    void compute();
};


                    /* TWO-DIMENSIONAL MODEL FUNCTION BASE CLASS */

class ModelFunction2D: public AbstractModelFunction
{
public:
    ModelFunction2D(user_func2D_t user_func);
    ModelFunction2D(user_func2D_t user_func, std::vector<double> &pars, void *extra_pars = nullptr);

    virtual ~ModelFunction2D();

    void setArgument(std::vector<double> &argX, std::vector<double> &argY);
    void setArgument(double xmin, double xmax, double ymin, double ymax, double xstep = 1.0, double ystep = 1.0);

    std::vector<double> operator()(std::vector<double> &x, std::vector<double> &y);

protected:
    user_func2D_t userFunction;
    std::vector<double> argumentX, argumentY;

    void compute();
};



                    /* GAUSSIAN AND MOFFAT MODEL FUNCTION CLASSES */

class GaussModelFunction: public ModelFunction
{
public:
    GaussModelFunction();
    GaussModelFunction(std::vector<double> &pars, void *extra_pars); // if extra_pars != nullptr then it is interpreted as a pointer to size_t
                                                                     // and *extra_pars is degree of background polynom

private:
    void checkParams();
    void checkConstrains();
};


class MoffatModelFunction: public ModelFunction
{
public:
    MoffatModelFunction();
    MoffatModelFunction(std::vector<double> &pars, void *extra_pars); // if extra_pars != nullptr then it is interpreted as a pointer to size_t
                                                                     // and *extra_pars is degree of background polynom

private:
    void checkParams();
    void checkConstrains();
};


class Gauss2DModelFunction: public ModelFunction2D
{
public:
    Gauss2DModelFunction();
    Gauss2DModelFunction(std::vector<double> &pars, void *extra_pars = nullptr); // if extra_pars != nullptr then it is interpreted as a pointer to size_t:
                                                                       // extra_pars[0] and extra_pars[1] are degrees of background polynom
                                                                       // along X and Y axis

private:
    void checkParams();
    void checkConstrains();
};


class Moffat2DModelFunction: public ModelFunction2D
{
public:
    Moffat2DModelFunction();
    Moffat2DModelFunction(std::vector<double> &pars, void *extra_pars = nullptr); // if extra_pars != nullptr then it is interpreted as a pointer to size_t:
                                                                        // extra_pars[0] and extra_pars[1] are degrees of background polynom
                                                                        // along X and Y axis

private:
    void checkParams();
    void checkConstrains();
};

#endif // MODELFUNCTION_H
