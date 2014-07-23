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

#ifndef CASADI_LAPACK_LU_DENSE_HPP
#define CASADI_LAPACK_LU_DENSE_HPP

#include "casadi/core/function/linear_solver_internal.hpp"
#include <casadi/interfaces/lapack/casadi_linearsolver_lapacklu_export.h>

namespace casadi {

/** \defgroup plugin_LinearSolver_lapacklu
*
   * This class solves the linear system <tt>A.x=b</tt> by making an LU factorization of A: \n
   * <tt>A = L.U</tt>, with L lower and U upper triangular
   *
*/ 

/** \pluginsection{LinearSolver,lapacklu} */

/// \cond INTERNAL

  /// LU-Factorize dense matrix (lapack)
  extern "C" void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

  /// Solve a system of equation using an LU-factorized matrix (lapack)
  extern "C" void dgetrs_(char* trans, int *n, int *nrhs, double *a,
                          int *lda, int *ipiv, double *b, int *ldb, int *info);

  /// Calculate col and row scaling
  extern "C" void dgeequ_(int *m, int *n, double *a, int *lda, double *r, double *c,
                          double *colcnd, double *rowcnd, double *amax, int *info);

  /// Equilibrate the system
  extern "C" void dlaqge_(int *m, int *n, double *a, int *lda, double *r, double *c,
                          double *colcnd, double *rowcnd, double *amax, char *equed);

  /** \brief \pluginbrief{LinearSolver,lapacklu}
   *
   * @copydoc LinearSolver_doc
   * @copydoc plugin_LinearSolver_lapacklu
   *
   */
  class CASADI_LINEARSOLVER_LAPACKLU_EXPORT LapackLUDenseInternal : public LinearSolverInternal {
  public:
    // Create a linear solver given a sparsity pattern and a number of right hand sides
    LapackLUDenseInternal(const Sparsity& sparsity, int nrhs);

    /** \brief  Create a new LinearSolver */
    static LinearSolverInternal* creator(const Sparsity& sp, int nrhs)
    { return new LapackLUDenseInternal(sp, nrhs);}

    /// Clone
    virtual LapackLUDenseInternal* clone() const;

    /// Destructor
    virtual ~LapackLUDenseInternal();

    /// Initialize the solver
    virtual void init();

    /// Prepare the solution of the linear system
    virtual void prepare();

    /// Solve the system of equations
    virtual void solve(double* x, int nrhs, bool transpose);

    /// A documentation string
    static const std::string meta_doc;

  protected:

    /// Scale columns
    void colScaling(double* x, int nrhs);

    /// Scale rows
    void rowScaling(double* x, int nrhs);

    // Matrix
    std::vector<double> mat_;

    /// Pivoting elements
    std::vector<int> ipiv_;

    /// Col and row scaling
    std::vector<double> r_, c_;

    /// Type of scaling during the last equilibration
    char equed_;

    /// Equilibrate?
    bool equilibriate_;

    /// Allow the equilibration to fail
    bool allow_equilibration_failure_;

    /// Dimensions
    int ncol_, nrow_;

  };

/// \endcond

} // namespace casadi

#endif // CASADI_LAPACK_LU_DENSE_HPP
