// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file NonlinearFunction.h
 * \brief Declare abstract base class NonlinearFunction.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURNONLINEARFUNCTION_H
#define MINOTAURNONLINEARFUNCTION_H

#include "Types.h"

namespace Minotaur {

  struct LTHessStor;
  class NonlinearFunction;
  class VarBoundMod;
  typedef boost::shared_ptr<NonlinearFunction> NonlinearFunctionPtr;
  typedef boost::shared_ptr<VarBoundMod> VarBoundModPtr;
  typedef std::vector<VarBoundModPtr> VarBoundModVector;
  typedef VarBoundModVector::iterator VarBoundModIter;

  /**
   * \brief Base class for nonlinear functions. 
   */
  class NonlinearFunction {

  public:
    /// Default constructor
    NonlinearFunction();

    /// Destroy.
    virtual ~NonlinearFunction();

    virtual void addConst(const double eps, int &err); 

    virtual void sqrRoot(int &err); 

    /** 
     * \brief Make a clone using new variables. 
     *
     * \param [in] vbeg it points to the variable id 0.  vbeg+k points to
     * variable id k, where k>=0. 
     * \param [out] err must be nonzero if function wasn't cloned.
     */
    virtual NonlinearFunctionPtr cloneWithVars(VariableConstIterator vbeg,
                                               int *err) const = 0;

    /**
     * \brief Calculate upper and lower bounds on the function using bounds of
     * the variables.
     * \param [out] lb Pointer that will contain the value of lower bound.
     * \param [out] lb Pointer that will contain the value of upper bound.
     * \param [out] error 0 if no error is encountered, nonzero otherwise.
     */
    virtual void computeBounds(double *lb, double *ub, int *error);

    /**
     * \brief Evaluate the function at a given point x.
     *
     * \param [in] x The size of the array x must exceed the highest index of the
     * variables used in the function.
     * \param [out] *error It should be set a positive value if there is
     * error encountered while evaluating. Leave undisturbed otherwise.
     * \return The value of function of x.
     */
    virtual double eval(const double *x, int *error) = 0;

    /**
     * \brief Evaluate and add gradient at a given point.
     *
     * \param [in] x The size of the array x must exceed the highest index of the
     * variables used in the function.
     * \param [out] grad_f The values of grad_f are incremented with the
     * gradient of this function at x. The array grad_f is dense.
     * \param [out] error It should be set a positive value if there is
     * error encountered while evaluating. Leave undisturbed otherwise.
     */
    virtual void evalGradient(const double *x, double *grad_f, int *error) 
      = 0;

    /**
     * \brief Evaluate and add hessian at a given point.
     *
     * \param [in] mult Multiplier for this objective/constraint function
     * \param [in] x The point where we need the hessian.
     * \param [in] stor The Hessian storage information.
     * \param [out] values The Hessian values to which we add Hessian of
     * this function.
     * \param [out] error We set it to nonzero if any errors are encountered. Do
     * not change it otherwise.
     */
    virtual void evalHessian(const double mult, const double *x, 
                             const LTHessStor *stor, double *values, 
                             int *error) = 0;

    /**
     * \brief Fill sparsity of hessian into hessian storage.
     *
     * \param [in/out] stor We add variables into stor->cols
     */
    virtual void  fillHessStor(LTHessStor *stor ) = 0;

    /**
     * \brief Evaluate and add gradient at a given point to the jacobian.
     *
     * \param [in] x The size of the array x must exceed the highest index of the
     * variables used in the function.
     * \param [out] values The values of jacobian are incremented with the
     * gradient of this function at x. 'values' only contains nonzeros of
     * jacobian. The indices (or offsets) where this nonlinear function
     * should put in the values should be calculated in the prepJac
     * function.  
     * \param [out] error It should be set a zero value if there is
     * error encountered while evaluating. Leave undisturbed otherwise.
     */
    virtual void fillJac(const double *x, double *values, int *error) = 0;

    /**
     * \brief Finalize hessian preparation. 
     *
     * \param [in] stor contains the sparsity pattern of hessian of
     * lagrangian. The nonlinear function should save offsets or make
     * other preparation to evaluate hessian.
     */
    virtual void  finalHessStor(const LTHessStor *stor) = 0;

    /**
     * \brief If a variable is fixed at a given value and removed, what is
     * the constant (offset) needed to be added.
     *
     * \param [in] v The variable that is fixed.
     * \param [in] val The value at which v is to be fixed.
     * \return Return the value of the offset.
     */
    virtual double getFixVarOffset(VariablePtr /* v */, double /* val */) 
    {assert (!"implment me!"); return 0;};

    /// Return the type of function: polynomial, ... 
    virtual FunctionType getType() const;

    /**
     * \brief Get variables used in this function
     *
     * \param [in] vars A set of variable-pointers into which variables are
     * to be inserted.
     */
    virtual void getVars(VariableSet *vars) = 0;

    /**
     * \brief Check if function contains a variable.
     *
     * \param [in] v The variable that we want to test.
     * \return True if this function is has v. False if it doesn't use it. 
     */
    virtual bool hasVar(ConstVariablePtr v) const;

    /**
     * \brief Check if the function is a sum of square terms
     *
     * \return True if this function is a sum of squares, False otherwise
     */
    virtual bool isSumOfSquares() const {
      return true;
    }

    /**
     * \brief Multiply by a constant.
     *
     * \param [in] c double value with which we want to multiply.
     */
    virtual void multiply(double c) = 0;

    /// Return the number of variables in this function.
    virtual UInt numVars() { return vars_.size(); };

    /**
     * \brief Prepare for evaluating sparse jacobian. 
     *
     * All the variables that are in this function can occur in one or more
     * of linear, quadratic and nonlinear functions. All variables that
     * occur in the whole function can be accessed by iterating between vbeg
     * and vend. This function is used to find the offset for variables that
     * occur in the nonlinear function.
     */
    virtual void prepJac(VarSetConstIter vbeg, VarSetConstIter vend) = 0;

    /// Remove a variable v from the function after fixing it to value val.
    virtual void removeVar(VariablePtr /* v */, double /* val */)
    {assert (!"implement me!");};

    /// Substitute a variable with another.
    virtual void subst(VariablePtr /* out */, VariablePtr /* in */,
                       double /* rat */)
    {assert (!"implement me!");};

    /**
     * \brief Take perspective of this function with respect to a given variable.
     *
     * Perspective of a given function f(x) with respect to a given variable z
     * is g(x,z) = z.f(x/z)
     * \param  [in] z The variable for which you take the perspective
     * \param  [in] eps The tolerance to tackle function value at 0 value. If
     * the perspective is approximated in the computational graph, the epsilon
     * is added to the z variable. If the perspective uses exact graph, then
     * an approximation in Hessian is used when z < eps.
     * \param [out] err must be nonzero if function wasn't cloned.
     * \return A new nonlinear function with an additional variable that gives
     * the perspective of this function 
     */
    virtual NonlinearFunctionPtr getPersp(VariablePtr z, double eps,
                                          int *err) const;

    virtual void varBoundMods(double /* lb */, double /* ub */,
                              VarBoundModVector & /* mods */,
                              SolveStatus * /* status */)
    {};

    /// \return first iterator for the variables in this function.
    virtual VariableSet::iterator varsBegin()
    {return vars_.begin();};

    /// \return last iterator for the variables in this function.
    virtual VariableSet::iterator varsEnd()
    {return vars_.end();};

    /// Display the nonlinear function.
    virtual void write(std::ostream &out) const;

  protected:
    /// A set of variables used in this function.
    VariableSet vars_;
  };

  typedef boost::shared_ptr<NonlinearFunction> NonlinearFunctionPtr;
  typedef boost::shared_ptr<const NonlinearFunction> ConstNonlinearFunctionPtr;  

}
#endif

// Local Variables: 
// mode: c++ 
// eval: (c-set-style "k&r") 
// eval: (c-set-offset 'innamespace 0) 
// eval: (setq c-basic-offset 2) 
// eval: (setq fill-column 78) 
// eval: (auto-fill-mode 1) 
// eval: (setq column-number-mode 1) 
// eval: (setq indent-tabs-mode nil) 
// End:
