// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file QuadraticFunction.h
 * \brief Construct and manage quadratic functions.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURQUADRATICFUNCTION_H
#define MINOTAURQUADRATICFUNCTION_H

#include "Types.h"

namespace Minotaur {
  class LinearFunction;
  struct LTHessStor;
  class QuadraticFunction;
  class PolynomialFunction;
  class Variable;
  typedef boost::shared_ptr<QuadraticFunction> QuadraticFunctionPtr;
  typedef boost::shared_ptr<const QuadraticFunction> ConstQuadraticFunctionPtr;
  typedef boost::shared_ptr<PolynomialFunction> PolyFunPtr;
  typedef boost::shared_ptr<const LinearFunction> ConstLinearFunctionPtr;
  typedef boost::shared_ptr<LinearFunction> LinearFunctionPtr;
  typedef boost::shared_ptr<Variable> VariablePtr;

  //class EigenVector;
  //typedef boost::shared_ptr<const EigenVector> EigenVectorPtr;
  //typedef std::pair<double, EigenVectorPtr> EigenPair;

  class QuadraticFunction {
    public:
      /// Default constructor
      QuadraticFunction();

      /// Constructor for a matrix.
      QuadraticFunction(UInt nz, double *vals, UInt *irow, UInt *jcol,
                        VariableConstIterator vbeg );

      /// Destroy
      ~QuadraticFunction();

      /** 
       * Copy the linear function. Variables and weights are copied. The weights
       * in the clone and the original do not share the same space in memory.
       */
      QuadraticFunctionPtr clone() const;

      QuadraticFunctionPtr cloneWithVars(VariableConstIterator vbeg) const;

      /**
       * Add a term of the form a*x_i*x_j to the expression. We do not check
       * if a similar term does not exist. 
       */
      void addTerm(VariablePair vp, const double weight); 

      /**
       * Add a term of the form a*x_i*x_j to the expression. We do not check
       * if a similar term does not exist. 
       */
      void addTerm(ConstVariablePtr v1, ConstVariablePtr v2, double weight);

      double getFixVarOffset(VariablePtr v, double val);

      /**
       * Add a term of the form a*x_i*x_j to the expression. If a similar term
       * exists, it is incremented. If the new coefficient becomes zero, the
       * term is dropped.
       */
      void incTerm(ConstVariablePair vp, const double weight);

      /**
       * Add a term of the form a*x_i*x_j to the expression. If a similar term
       * exists, it is incremented. If the new coefficient becomes zero, the
       * term is dropped.
       */
      void incTerm(ConstVariablePtr v1, ConstVariablePtr v2, const double weight);

      /*
       * Multiply by a constant. If constant is zero, all terms are removed.
       */
      void mult(double c);

      /**
       * Remove a variable v from the function. Add to lf any linear terms
       * obtained by fixing the variable to the value 'val'.
       */
      void removeVar(VariablePtr v, double val, LinearFunctionPtr lf);

      void subst(VariablePtr out, VariablePtr in, double rat);

      /**
       * Add a quadratic function to this quadratic function. Terms that become 
       * zero are not retained in the function.
       */
      friend QuadraticFunctionPtr operator + (ConstQuadraticFunctionPtr q1, 
          ConstQuadraticFunctionPtr q2);

      /**
       * Subtract a linear function from this function. Terms that become zero
       * are still retained in the function.
       */
      friend QuadraticFunctionPtr operator-(ConstQuadraticFunctionPtr q1, 
          ConstQuadraticFunctionPtr q2);

      /// Multiply a quadratic function with a constant.
      friend QuadraticFunctionPtr operator*(const double c, 
          ConstQuadraticFunctionPtr q2);

      /// Multiply a linear function and quadratic function.
      friend PolyFunPtr operator*(ConstQuadraticFunctionPtr q2, 
          ConstLinearFunctionPtr l1);

      /// Multiply two quadratics
      friend PolyFunPtr operator*(ConstQuadraticFunctionPtr q1, 
          ConstQuadraticFunctionPtr q2);


      /**
       * This increment operator is dangerous to use because it only works on
       * objects of type QuadraticFunction and does not work on type
       * QuadraticFunctionPtr.  So if you have:
       * 
       * QuadraticFunctionPtr qPtr, q2Ptr;
       * 
       * Then you cannot do:
       * 
       * qPtr += q2Ptr;
       * 
       * You will have to do:
       * 
       * (*qPtr) += q2Ptr;
       * 
       * The user must check if the left operand is not NULL.
       */
      void operator+=(ConstQuadraticFunctionPtr q2);

      /**
       * Multiply the quadratic with a constant. Same precaution as for +=
       * operator above. If c is zero, then quadratic becomes empty. It is
       * better for the calling routine to check if c is zero, if so, just
       * delete the quadratic.
       */
      void operator*=(const double c);

      /// Get the list of variables and how many times they occur
      VarCountConstMap * getVarMap() const;

      /// Constant Iterators to visit all the quadratic terms: begin.
      VariablePairGroupConstIterator begin() const { return(terms_.begin()); }

      /// Constant Iterators to visit all the quadratic terms: end.
      VariablePairGroupConstIterator end() const { return(terms_.end()); }

      /**
       * Get the coefficient of a term. Returns 0 if the term does not exist
       * in the function.
       */
      double getWeight(ConstVariablePair & vp) ;

      /**
       * Get the coefficient of a term. Returns 0 if the term does not exist
       * in the function.
       */
      double getWeight(ConstVariablePtr v1, ConstVariablePtr v2) ;

      /// Get the number of times variable v1 occurs in quadratic terms.
      int getFreq(ConstVariablePtr v1);

      /**
       * \brief Check if function contains a variable.
       *
       * \param[in] v The variable that we want to test.
       * \return True if this function is has v. False if it doesn't use it. 
       */
      bool hasVar(ConstVariablePtr v) const;

      /**
       * Evaluate the value of this quadratic expression for a given point x.
       * It is assumed that x[i] is the value of the variable whose
       * var->getIndex() returns i.
       */
      double eval(const std::vector<double> &x) const;

      /**
       * Evaluate the value of this quadratic expression for a given point x.
       * It is assumed that x[i] is the value of the variable whose
       * var->getIndex() returns i.
       */
      double eval(const double *x) const;

      /**
       * Evaluate the values of the gradient of the quadratic expression at a
       * given point x.  It is assumed that x[i] is the value of the variable
       * whose var->getIndex() returns i.
       */
      void evalGradient(const double *x, double *grad_f); 

      /**
       * Evaluate the values of the gradient of the quadratic expression at a
       * given point x.  It is assumed that x[i] is the value of the variable
       * whose var->getIndex() returns i.
       */
      void evalGradient(const std::vector<double> &x,
                        std::vector<double> &grad_f); 
      void evalHessian(const double mult, const double *x, 
                       const LTHessStor *stor, double *values , int *error);

      void prepJac(VarSetConstIter vbeg, VarSetConstIter vend);
      void prepHess();

      void fillHessStor(LTHessStor *hess);
      void fillJac(const double *x, double *values, int *error);
      void finalHessStor(const LTHessStor *hess);

      /// Get the number of terms in this expression
      UInt getNumTerms() const;

      /// Get the number of variables in this expression
      UInt getNumVars() const;

      void getVars(VariableSet *vars);

      /**
       * Return true if the quadratic expression always returns zero on
       * evaluation.
       */
      bool isZero() const { return (getNumTerms() == 0); }

      /// Iterate over the variables in the quadratic function: begin.
      VarIntMapConstIterator varsBegin() const;

      /// Iterate over the variables in the quadratic function: end.
      VarIntMapConstIterator varsEnd() const;

      /// Display the quadratic function
      void write(std::ostream &s) const;

    private:
      /// Tolerance below which a coefficient is deemed zero
      const double etol_;

      double *hCoeffs_;
      UInt *hFirst_;
      UInt *hOff_;
      UInt *hSecond_;
      UIntVector hStarts_;

      UIntVector jacOff_;
      UIntVector jacInd_;

      /**
       * Each item in this map is a pair of a variable-pair and a
       * weight. A variable pair is in turn a pair of two variables.
       */
      VariablePairGroup terms_;

      /**
       * Set of variables that are in this quadratic and in how many times
       * do they occur. e.g. x0 occurs 4 times in:
       * x0^2 + 2x0x1 + x0x3
       */
      VarIntMap varFreq_;

      void sortLT_(UInt n, UInt *f, UInt *s, double *c);
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
