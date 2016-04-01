// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file LinearFunction.h
 * \author Ashutosh Mahajan, Argonne National Laboratory.
 * \brief Declare the class LinearFunction for storing and modifying a 
 * linear function.
 */

#ifndef MINOTAURLINEARFUNCTION_H
#define MINOTAURLINEARFUNCTION_H

#include "Types.h"

namespace Minotaur {
  class LinearFunction;
  class QuadraticFunction;
  class Variable;

  typedef boost::shared_ptr<LinearFunction> LinearFunctionPtr;
  typedef boost::shared_ptr<const LinearFunction> ConstLinearFunctionPtr;
  typedef boost::shared_ptr<QuadraticFunction> QuadraticFunctionPtr;
  typedef boost::shared_ptr<const QuadraticFunction> ConstQuadraticFunctionPtr;
  typedef boost::shared_ptr<Variable> VariablePtr;


  /// The base class linear function is of the form c'x.
  class LinearFunction {
  public:
    /// Default constructor
    LinearFunction();

    /** 
     * Constructor with a tolerance level below which a coefficient is
     * considered to be zero.
     */
    LinearFunction(const double tol);

    /**
     * Construct a linear function from a coefficient array of size n and n
     * variables.
     */
    LinearFunction(double *a, VariableConstIterator vbeg, 
                   VariableConstIterator vend, double tol); 


    /// Destroy
    ~LinearFunction();

    void add(LinearFunctionPtr lf);

    /**
     * Add new a linear term to this linear function, with coefficient a. Use
     * this method only when you are sure that the linear function does not
     * already contain this variable. Otherwise use incTerm(). If the weight "a" 
     * is zero, then nothing is added.
     */
    void addTerm(ConstVariablePtr var, const double a); 

    /**
     * Removes all terms from the function
     */
    void clearAll();

    /** 
     * Copy the linear function. Variables and weights are copied. The weights
     * in the clone and the original do not share the same space in memory.
     */
    LinearFunctionPtr clone() const;

    LinearFunctionPtr cloneWithVars(VariableConstIterator vbeg) const;

    /**
     * \brief Get bounds based on lower and upperbounds of each variable
     * \param [out] l This pointer should contain lower bound.
     * \param [out] u This pointer should contain upper bound.
     */
    void computeBounds(double *l, double *u);

    /**
     * Evaluate the value of this linear function at a given vector of
     * doubles. It is assumed that the x[i] has the value for variable with id
     * 'i'.
     */
    double eval(const std::vector<double> &x) const;

    /**
     * Evaluate the value of this linear function at a given vector of
     * doubles. It is assumed that the x[i] has the value for variable with id
     * 'i'.
     */
    double eval(const double *x) const;

    /**
     * Evaluate the gradient of this linear function. It is assumed that x[id]
     * will have the gradient along the direction of variable with ID=id.
     */
    void evalGradient(double *grad_f) const;

    void fillJac(double *values, int *error);

    double getFixVarOffset(VariablePtr v, double val);

    /// Get the number of terms in this function.
    UInt getNumTerms() const { return(terms_.size()); }

    void getVars(VariableSet *vars);

    /**
     * Get the weight of a variable in this function. If the variable is not
     * found, it returns zero. Conversely, if the weight returned is zero,
     * then variable is not stored in the data structures of this function. 
     */
    double getWeight(ConstVariablePtr var) const;

    /**
     * \brief Check if function contains a variable.
     *
     * \param[in] v The variable that we want to test.
     * \return True if this function is has v. False if it doesn't use it. 
     */
    bool hasVar(ConstVariablePtr v) const;

    /**
     * Add new a linear term to this linear function, with coefficient a. If
     * the function already contains this variable, then the value is
     * incremented. If the new value becomes zero, the variable is dropped.
     */
    void incTerm(ConstVariablePtr var, const double a);

    /// Multiply the linear function by a number.
    void multiply(double d);

    void prepJac(VarSetConstIter vbeg, VarSetConstIter vend);

    /// Remove a variable v from the function.
    void removeVar(VariablePtr v, double val);

    /// Square the function.
    QuadraticFunctionPtr square();

    /// Iterate over the terms in the linear function: begin.
    VariableGroupConstIterator termsBegin() const;

    /// Iterate over the terms in the linear function: end.
    VariableGroupConstIterator termsEnd() const;

    /// Writes the function to a stream.
    void write(std::ostream &out) const;


    /**
     * Add a linear function to this linear function. Terms that become zero
     * are still retained in the function.
     */
    friend LinearFunctionPtr operator + (ConstLinearFunctionPtr l1, 
                                         ConstLinearFunctionPtr l2);

    /**
     * Subtract a linear function from this function. Terms that become zero
     * are still retained in the function.
     */
    friend LinearFunctionPtr operator-(ConstLinearFunctionPtr l1, 
                                       ConstLinearFunctionPtr l2);

    /// Multiply a linear function with a constant.
    friend LinearFunctionPtr operator*(const double c, 
                                       ConstLinearFunctionPtr l2);

    /**
     * Multiply two linear functions to get a quadratic function. If either of
     * the linear functions is NULL, return NULL.
     */
    friend QuadraticFunctionPtr operator*(ConstLinearFunctionPtr l1, 
                                          ConstLinearFunctionPtr l2);

    /**
     * This increment operator is dangerous to use because it only works on
     * objects of type LinearFunction and does not work on type
     * LinearFunctionPtr.  So if you have:
     * 
     * LinearFunctionPtr lPtr, l2Ptr;
     * 
     * Then you cannot do:
     * 
     * lPtr += l2Ptr;
     * 
     * You will have to do:
     * 
     * (*lPtr) += l2Ptr;
     * 
     * The user must ensure left operand is not NULL.
     */
    void operator+=(ConstLinearFunctionPtr l2);

    /** 
     * Subtract l2 from this linear function.
     * The user mu ensure left operand is not NULL.
     */
    void operator-=(ConstLinearFunctionPtr l2);

    /**
     * Multiply with a constant. Same precaution as for +=
     * operator above. If c is zero, then function becomes empty. It is
     * better for the calling routine to check if c is zero, if so, just
     * delete the function.
     */
    void operator*=(const double c);


  private:
    /**
     * True if terms in linear function are modified since previous call to
     * prepJac.
     */
    bool hasChanged_;

    /// Offsets for jacobian.
    DoubleVector off_;

    /**
     * terms_ is a map with variables as keys and their coefficients as
     * values.
     */
    VariableGroup terms_;

    /// Tolerance below which a coefficient is considered 0.
    double tol_;

    /// Copy constructor is not allowed.
    LinearFunction (const LinearFunction &l);

    /// Copy by assignment is not allowed.
    LinearFunction  & operator = (const LinearFunction &l);

  };
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
