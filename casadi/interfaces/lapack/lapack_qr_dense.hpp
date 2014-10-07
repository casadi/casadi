/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef CASADI_LAPACK_QR_DENSE_HPP
#define CASADI_LAPACK_QR_DENSE_HPP

#include "casadi/core/function/linear_solver_internal.hpp"
#include <casadi/interfaces/lapack/casadi_linearsolver_lapackqr_export.h>

/** \defgroup plugin_LinearSolver_lapackqr
*
* This class solves the linear system <tt>A.x=b</tt> by making an QR factorization of A: \n
* <tt>A = Q.R</tt>, with Q orthogonal and R upper triangular
*/

/** \pluginsection{LinearSolver,lapackqr} */

/// \cond INTERNAL
namespace casadi {

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

  /** \brief  \pluginbrief{LinearSolver,lapackqr}
   *
   @copydoc LinearSolver_doc
   @copydoc plugin_LinearSolver_lapackqr
   *
   */
  class CASADI_LINEARSOLVER_LAPACKQR_EXPORT LapackQrDense : public LinearSolverInternal {
  public:
    // Create a linear solver given a sparsity pattern and a number of right hand sides
    LapackQrDense(const Sparsity& sparsity, int nrhs);

    // Clone
    virtual LapackQrDense* clone() const;

    /** \brief  Create a new LinearSolver */
    static LinearSolverInternal* creator(const Sparsity& sp, int nrhs)
    { return new LapackQrDense(sp, nrhs);}

    // Destructor
    virtual ~LapackQrDense();

    // Initialize the solver
    virtual void init();

    // Prepare the solution of the linear system
    virtual void prepare();

    // Solve the system of equations
    virtual void solve(double* x, int nrhs, bool transpose);

    /// A documentation string
    static const std::string meta_doc;

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

} // namespace casadi

/// \endcond
#endif // CASADI_LAPACK_QR_DENSE_HPP
