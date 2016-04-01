// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file PolynomialFunction.h
 * \brief Define the MonomialFunction and PolynomialFunction classs for
 * handling polynomial constraints.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURPOLYNOMIALFUNCTION_H
#define MINOTAURPOLYNOMIALFUNCTION_H

#include "NonlinearFunction.h"

namespace Minotaur {

  class CGraph;
  class CNode;
  class LinearFunction;
  class MonomialFunction;
  class QuadraticFunction;
  typedef boost::shared_ptr<CGraph> CGraphPtr;
  typedef boost::shared_ptr<LinearFunction> LinearFunctionPtr;
  typedef boost::shared_ptr<MonomialFunction> MonomialFunPtr;
  typedef boost::shared_ptr<QuadraticFunction> QuadraticFunctionPtr;
  typedef boost::shared_ptr<const LinearFunction> ConstLinearFunctionPtr;
  typedef boost::shared_ptr<const MonomialFunction> ConstMonomialFunPtr;
  typedef boost::shared_ptr<const QuadraticFunction> ConstQuadraticFunctionPtr;
  typedef std::vector<MonomialFunPtr> MonomialVector;
  typedef MonomialVector::iterator MonomialIter;
  typedef MonomialVector::const_iterator MonomialConstIter;


  class PolynomialFunction;
  typedef boost::shared_ptr<PolynomialFunction> PolyFunPtr;
  typedef boost::shared_ptr<const PolynomialFunction> ConstPolyFunPtr;  

  /**
   * \brief MonomialFunction represents functions of the form
   * of the form \f$\prod_ia_ix_i^{b_i}\f$ where \f$a_i \in \mathbb{R}, b_i
   * \in \mathbb{N}\f$.
   */
  class MonomialFunction : public NonlinearFunction {
  public:
    /// Default constructor.
    MonomialFunction();

    /// A constant monomial.
    MonomialFunction(double c);

    /// Constructor with a single variable, a  multiplier, and a power.
    MonomialFunction(double c, ConstVariablePtr v, UInt p);

    /// Destroy.
    ~MonomialFunction();

    /// Create a full copy.
    MonomialFunPtr clone() const;

    /// Make a clone using new variables. vbeg points to the variable id 0.
    /// vbeg+k points to variable id k, where k>=0.
    NonlinearFunctionPtr cloneWithVars(VariableConstIterator vbeg, 
                                       int *err) const;

    // Evaluate.
    double eval(const double *x, int *error);

    // Evaluate.
    double eval(const DoubleVector &x, int *error);

    // Gradient.
    void evalGradient(const double *x, double *grad_f, int *error);

    void evalHessian(const double , const double *, 
                     const LTHessStor *, double *, 
                     int *) 
    {assert(!"can't fill hessian in a monomial function.");};

    /// Populate the cgraph with a subgraph representing this monomial and
    /// return the root node of this subgraph.
    CNode* fillCG(CGraphPtr cg);

    void fillJac(const double *, double *, int *) 
    {assert(!"can't fill jacobian in a monomial function.");};

    void  fillHessStor(LTHessStor *)
    {assert(!"can't fill hessian in a monomial function.");};

    void  finalHessStor(const LTHessStor *) 
    {assert(!"can't fill hessian in a monomial function.");};

    double getCoeff() const;

    /// Get the degree, sum of all powers.
    int getDegree() const;

    /// Monomial terms
    const VarIntMap* getTerms() const;

    void getVars(VariableSet *) {};

    /// Multiply with a variable raised to power.
    void multiply(double coeff, ConstVariablePtr v, int p);

    /// Multiply with another monomial.
    void multiply(ConstMonomialFunPtr m2);

    /// Multiply a constant.
    void multiply(double c);

    void prepJac(VarSetConstIter, VarSetConstIter)
    {assert(!"can't fill jacobian in a monomial function.");};

    /// Iterator for the first term.
    VarIntMapConstIterator termsBegin();

    /// Iterator at the end of the terms.
    VarIntMapConstIterator termsEnd();

    /// Raise the monomial to power k.
    void toPower(int k);

    /// Display.
    void write(std::ostream &out) const;

    /// Multiply a constant.
    void operator*=(double c);

    /// Multiply a monomial.
    void operator*=(ConstMonomialFunPtr m2);

  private:
    /// The coefficient.
    double coeff_;

    /// The degree of monomial.
    int deg_;

    /**
     * \brief If a coefficient's magnitude is lower than eTol_, it is considered
     * zero.
     */
    const double eTol_;

    /// The variables and their exponents.
    VarIntMap terms_;
  };

  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------

  /**
   * \brief PolynomialFunction represents functions of the form
   * \f$\sum_jm_j\f$, where \f$m_j\f$ is a MonomialFunction.
   */
  class PolynomialFunction : public NonlinearFunction {
  public:
    /// Default constructor.
    PolynomialFunction();

    /// Expand polynomial from a computational graph.
    PolynomialFunction(CGraphPtr cg);

    /// Destroy
    ~PolynomialFunction();

    /// Make a clone.
    PolyFunPtr clone() const;

    /// Add a monomial m to this polynomial. The monomial is first cloned.
    void add(ConstMonomialFunPtr m);

    // base class function.
    NonlinearFunctionPtr cloneWithVars(VariableConstIterator vbeg, 
                                       int *err) const;

    /// Create a cgraph representing this polynomial and save it as a private
    /// member. Useful for solving using native derivatives.
    void createCG();

    // base class function.
    double eval(const double *x, int *error);

    // base class function.
    void evalGradient(const double *x, double *grad_f, int *error);

    // base class function.
    void evalHessian(const double mult, const double *x, 
                     const LTHessStor *stor, double *values, 
                     int *error);

    // base class function.
    void  fillHessStor(LTHessStor *stor);

    // base class function.
    void fillJac(const double *x, double *values, int *error);

    // base class function.
    void  finalHessStor(const LTHessStor *stor);

    // base class function.
    void getVars(VariableSet *);

    /// Get the constant term.
    double getConstant();

    /// Return type of function.
    FunctionType getType() const;

    /// Return true if the constant is 0 and there are no terms.
    bool isEmpty() const;

    /// Multiply with (lf + c), where lf is a linear function.
    void multiply(ConstLinearFunctionPtr lf, double c);

    /// Multiply with c, a constant.
    void multiply(double c);

    // base class function.
    void prepJac(VarSetConstIter vb, VarSetConstIter ve);

    void recCG_(const CNode* cnode, double *c, MonomialVector *terms);
    void recCGMult_(MonomialVector *t1, MonomialVector *t2,
                    double c1, double c2,
                    MonomialVector *terms, double *c);
    /**
     * \brief If the polynomial has a constant term, return it, and remove it
     * from the polynomial.
     */
    double removeConstant();

    /**
     * If the polynomial has a linear function, add it to 'lf', and remove
     * it from the polynomial. If lf is empty and polynomial has linear
     * terms, new space is allocated.
     */
    void removeLinear(LinearFunctionPtr lf);

    /**
     * If the polynomial has a quadratic function, add it to 'qf', and remove
     * it from the polynomial. If qf is NULL and polynomial has quadratic
     * terms, new space is allocated.
     */
    void removeQuadratic(QuadraticFunctionPtr qf);

    /// Get an iterator to the first monomial
    MonomialConstIter termsBegin();

    /// Get an iterator to the last monomial
    MonomialConstIter termsEnd();

    /// Display.
    void write(std::ostream &out) const;

    /**
     * Add a monomial to this polynomial. The monomial is copied. The user
     * can change the monomial later wihout affecting the polynomial. This
     * function is dangerous because it operates on the raw pointer of the
     * polynomial and not the shared pointer. The user must ensure that it
     * is not null.
     */
    void operator+=(ConstMonomialFunPtr m);

    /// Add a polynomial to the polynomial.
    void operator+=(ConstPolyFunPtr p);

    /// Multiply by a constant.
    void operator*=(double c);

    /// Add a constant.
    void operator+=(double c);

    /// Add a linear function.
    void operator+=(ConstLinearFunctionPtr lf);

    /// Subtract a linear function.
    void operator-=(ConstLinearFunctionPtr lf);

    /// Multiply with a linear function.
    void operator*=(ConstLinearFunctionPtr lf);

    /// Add a quadratic.
    void operator+=(ConstQuadraticFunctionPtr qf);

    /// Multiply with a quadratic.
    void operator*=(ConstQuadraticFunctionPtr qf);

    /// Multiply with a polynomial.
    void operator*=(ConstPolyFunPtr p2);

    /// Add two polynomials. Creates a new polynomial.
    friend PolyFunPtr operator + (ConstPolyFunPtr p1, ConstPolyFunPtr p2);

    /// Subtract two polynomials.
    friend PolyFunPtr operator - (ConstPolyFunPtr p1, ConstPolyFunPtr p2);

    /// Multiply a polynomial with constant.
    friend PolyFunPtr operator*(double c, ConstPolyFunPtr p2);

    /// Multiply a linear function and polynomial function.
    friend PolyFunPtr operator*(ConstPolyFunPtr p2, 
                                ConstLinearFunctionPtr l1);

    /// Multiply a polynomial with quadratic.
    friend PolyFunPtr operator*(ConstPolyFunPtr p1, 
                                ConstQuadraticFunctionPtr q2);

    /// Multiply a polynomial with constant.
    friend PolyFunPtr operator*(ConstPolyFunPtr p1, ConstPolyFunPtr p2);

  protected:
    /// Additional constant.
    double cb_;

    /// If the polynomial is constructed from a CGraph, we keep a pointer.
    CGraphPtr cg_;

    /// Tolerance.
    const double eTol_;

    /// Each monomial term.
    MonomialVector terms_;

    /// Clear/reset all terms. Polynomial becomes 0. 
    void clear_();
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
