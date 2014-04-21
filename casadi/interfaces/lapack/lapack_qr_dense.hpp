/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef LAPACK_QR_DENSE_HPP
#define LAPACK_QR_DENSE_HPP

#include "casadi/symbolic/function/linear_solver_internal.hpp"
#include <casadi/interfaces/lapack/casadi_lapack_interface_export.h>

namespace casadi{

  /** \brief  Forward declaration of internal class
      @copydoc LinearSolver_doc
  */
  class LapackQRDenseInternal;

  /** \brief  QR LinearSolver with Lapack Interface
   *
   @copydoc LinearSolver_doc
   *
   * This class solves the linear system <tt>A.x=b</tt> by making an QR factorization of A: \n
   * <tt>A = Q.R</tt>, with Q orthogonal and R upper triangular
   *
   * LapackQRDense is an casadi::Function mapping from 2 inputs
   * [ A (matrix),b (vector)] to one output [x (vector)].
   *
   * The usual procedure to use LapackQRDense is: \n
   *  -# init()
   *  -# set the first input (A)
   *  -# prepare()
   *  -# set the second input (b)
   *  -# solve()
   *  -# Repeat steps 4 and 5 to work with other b vectors.
   *
   * The method evaluate() combines the prepare() and solve() step and is
   * therefore more expensive if A is invariant.
   *
   */
  class CASADI_LAPACK_INTERFACE_EXPORT LapackQRDense : public LinearSolver{
  public:

    /// Default (empty) constructor
    LapackQRDense();

    /// Create a linear solver given a sparsity pattern
    explicit LapackQRDense(const Sparsity& sparsity, int nrhs=1);

    /// Access functions of the node
    LapackQRDenseInternal* operator->();
    const LapackQRDenseInternal* operator->() const;

    /// Static creator function
#ifdef SWIG
    %callback("%s_cb");
#endif
    static LinearSolver creator(const Sparsity& sp, int nrhs){ return LapackQRDense(sp, nrhs);}
#ifdef SWIG
    %nocallback;
#endif

  };

/// \cond INTERNAL
#ifndef SWIG

  /// QR-factorize dense matrix (lapack)
  extern "C" void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau,
                          double *work, int *lwork, int *info);

  /// Multiply right hand side with Q-transpose (lapack)
  extern "C" void dormqr_(char *side, char *trans, int *n, int *m, int *k, double *a,
                          int *lda, double *tau, double *c, int *ldc,
                          double *work, int *lwork, int *info);

  /// Solve upper triangular system (lapack)
  extern "C" void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n,
                         double *alpha, double *a, int *lda, double *b, int *ldb);

  /// Internal class
  class CASADI_LAPACK_INTERFACE_EXPORT LapackQRDenseInternal : public LinearSolverInternal{
  public:
    // Create a linear solver given a sparsity pattern and a number of right hand sides
    LapackQRDenseInternal(const Sparsity& sparsity, int nrhs);

    // Clone
    virtual LapackQRDenseInternal* clone() const;

    // Destructor
    virtual ~LapackQRDenseInternal();

    // Initialize the solver
    virtual void init();

    // Prepare the solution of the linear system
    virtual void prepare();

    // Solve the system of equations
    virtual void solve(double* x, int nrhs, bool transpose);

  protected:

    // Matrix
    std::vector<double> mat_;

    // The scalar factors of the elementary reflectors
    std::vector<double> tau_;

    // qr work array
    std::vector<double> work_;

    // Dimensions
    int ncol_, nrow_;

  };

#endif // SWIG
/// \endcond

} // namespace casadi



#endif //LAPACK_QR_DENSE_HPP

