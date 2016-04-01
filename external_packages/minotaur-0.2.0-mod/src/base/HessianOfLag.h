//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file HessianOfLag.h
 * \brief Declare class HessianOfLag to get information about hessian of
 * Lagrangean of a Problem.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURHESSIAN_H
#define MINOTAURHESSIAN_H

#include "Types.h"

namespace Minotaur {
  class Function;
  typedef boost::shared_ptr<Function> FunctionPtr;

  struct LTHessStor {
    UInt nz;
    UInt nlVars;
    VariablePtr *rows;
    std::deque<UInt> *colQs;
    UInt *cols;
    UInt *starts;
  };
  

  /**
   * Given a function f(x), the hessian of f(x) is defined as the matrix whose
   * element \f$(i,j)\f$ is: d^2f/dx_idx_j], i=0,1,...,n-1, j=0,1,...,n-1.
   * Hessian is thus a square symmetric matrix. We store only the
   * lower-triangular part of the matrix, i.e., i=0,1,...,n-1, j=0,...,i
   * 
   * Thus each function can have an associated Hessian. NLP solvers do not
   * usually require hessian of the objective function or a function from a
   * single constraint. They usually need hessian of the lagrangean function.
   * This function is defined as:
   * \f[ l = \lambda_0f + \sum_{i=1}^m\lambda_ig_i \f]
   * 
   * This class provides methods for building and managing hessian of the
   * lagrangean by calling the methods of the Functions in objective and
   * constraints.
   */
  class HessianOfLag {
    public:
      /// Default constructor
      HessianOfLag();

      HessianOfLag(Problem *p);

      /**
       * Construct Hessian of Lagrangean from an objective, obj,
       * and a vector of constraints cons. 'n' is the number of variables and
       * is also the size of square matrix.
       */
      //HessianOfLag(ConstObjPtr obj, const ConstraintVector & cons, 
      //    const UInt n);

      /// Destroy.
      virtual ~HessianOfLag();

      /**
       * Return the number of non-zero elements in the lower-triangular part
       * of the hessian, including the diagonal).
       */
      virtual UInt getNumNz() const;

      /**
       * Given arrays iRow and jCol of size getNumNz(), fill the values of
       * rows and columns that have non-zeros.
       * 
       * e.g. if f(x) = x0^2 + x1^3x2 + x1 
       * 
       * then getNumNz() = 3
       *      iRow = [0, 1, 2]
       *      jCol = [0, 1, 1]
       */
      virtual void fillRowColIndices(UInt *irow, UInt *jcol);

      /**
       * Given arrays iRow and jCol of size getNumNz(), fill the values of
       * hessian matrix at a given point 'x', into the array 'values'. for the
       * above example, 
       * values = [2, 6x1x2, 3x1^2]
       */
      virtual void fillRowColValues(const double *x, double obj_mult,
                                    const double *con_mult, double *values, 
                                    int *error);

      /// Ugly hack to solve maximization problem. TODO: delete it.
      virtual void negateObj() {};

      virtual void setupRowCol();

      virtual void write(std::ostream &out) const;

    private:
      /** 
       * If a lagrange multiplier is less than etol_, then it is considered
       * zero.
       */
      double etol_;

      FunctionPtr obj_;

      /**
       * Stores for each constraint, where each hessian non-zero entry goes.
       * e.g.
       * x_1^2 + x_4x_5 <= 0
       * x_1^2 + x_2^2 + x_3^2 - 10*x_4x_5 <= 0
       * Then the non-zeros are stored for elements
       * (1,1), (2,2), (3,3), (4,5)
       * The offsets matrix is much bigger than the non-zeros. It contains:
       * 0, 3, 0, 1, 2, 3
       * The first 0 means that the hessian value of x^1 goes in position 0.
       * The hessian value for x_4x_5 goes in position 3 and so on.
       */

      Problem *p_;
      LTHessStor stor_;

  };

  typedef boost::shared_ptr<HessianOfLag> HessianOfLagPtr;
  typedef boost::shared_ptr<const HessianOfLag> ConstHessianOfLagPtr;  

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
