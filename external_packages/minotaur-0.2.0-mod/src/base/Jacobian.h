//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file Jacobian.h
 * \brief Get information about Jacobian of a given Problem.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURJACOBIAN_H
#define MINOTAURJACOBIAN_H

#include "Types.h"

namespace Minotaur {


  /**
   * This class is used for the Jacobian of a Problem. When a problem has
   * constraints of the form: 
   * \f[ g_i(x) \leq 0 \f] 
   * The Jacobian at a point
   * \f$x^*\f$ is a matrix whose \f$i-\f$th row is the gradient, 
   * \f$\nabla f_i(x^*)\f$.
   * 
   * This class also provides the sparsity pattern of the non-zeros in the
   * Jacobian. By default, we save the Jacobian constraint-wise and then
   * variable-wise. Thus, if we a have the constraints:
   * 
   * x_0^2 + x1^3 \leq 3
   * x_0 + x_2 + x_3 \leq 4
   * x_1 + x_3 \leq 4
   * x_3^3 \leq 8
   * 
   * We will save the indices of nonzeros in a 2-d matrix, rows corresponding 
   * to constraints and columns corresponding to variables. The matrix will
   * have entries:
   * 0 1
   * 0 2 3
   * 1 3
   * 3
   * 
   * The above scheme makes it easy to fill nonzeros in Constraint-Ordered
   * manner, but more expensive to fill in Variable-Ordered manner. So we also
   * use another 2-d matrix as above, but the matrix now looks like:
   * 0 0 1 1 2
   * 1 2   2
   *       3
   */
  class Jacobian {

    public:

      /// Default constructor.
      Jacobian();

      /// Constructor using constraints.
      Jacobian(const std::vector<ConstraintPtr> & cons, const UInt n);

      /// Destroy.
      virtual ~Jacobian();

      /// Return the number of nonzeros in the Jacobian.
      virtual UInt getNumNz();

      /**
       * Given arrays iRow and jCol, fill in the row and column index of each
       * non-zero in the jacobian.
       * e.g.
       * \f{eqnarray*} 
       * x_1x_3 + x_1^2x_4 &=& 2 \\ x_0 + x_5 &=& 0
       * \f}
       * 
       * then, iRow = [1 0 0 0 1], jCol = [0 1 3 4 5]. These indices are
       * arranged first in the order of increasing jCol and then increasing
       * iRow. 
       */
      virtual void fillRowColIndices(UInt *iRow, UInt *jCol); 

      /**
       * Given arrays iRow and jCol as above, fill in the row and column index
       * of each non-zero in the jacobian. For the above example, when \f$x =
       * (-1, 1, 1, 1, 0)\f$, then values = [3 -1 1 1 1].  
       */
      virtual void fillRowColValues(const double *x, double *values, 
          int *error);

      /// Fill indices, column wise.
      virtual void fillColRowIndices(UInt *, UInt *)
      { assert(!"implement me!");}
         
      /// Fill values, column wise.
      virtual void fillColRowValues(const double *, double *, int *)
      { assert(!"implement me!");}
         
      void write(std::ostream &out) const;

    private:
      /**
       * The constraints that constitute the system whose jacobian needs to be
       * evaluated.
       */
      const std::vector<ConstraintPtr> * cons_;

      /// Number of nonzeros
      UInt nz_;

  };
  typedef boost::shared_ptr<Jacobian> JacobianPtr;
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
