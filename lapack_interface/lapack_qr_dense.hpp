/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#include "casadi/fx/linear_solver.hpp"

namespace CasADi{
  
/** \brief  Forward declaration of internal class */
class LapackQRDenseInternal;

/** \brief  Public class */
class LapackQRDense : public LinearSolver{
public:

  /// Default (empty) constructor
  LapackQRDense();
  
  /// Create a linear solver given a sparsity pattern
  LapackQRDense(int nrow, int ncol, const std::vector<int>& rowind, const std::vector<int>& col, int nrhs=1);
    
  /// Access functions of the node
  LapackQRDenseInternal* operator->();
  const LapackQRDenseInternal* operator->() const;
};

/// QR-factorize dense matrix (lapack)
extern "C" void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

/// Multiply right hand side with Q-transpose (lapack)
extern "C" void dormqr_(char *side, char *trans, int *n, int *m, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *lwork, int *info);

/// Solve upper triangular system (lapack)
extern "C" void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);

/// Internal class
class LapackQRDenseInternal : public LinearSolverInternal{
  public:
    // Create a linear solver given a sparsity pattern and a number of right hand sides
    LapackQRDenseInternal(int nrow, int ncol, const std::vector<int>& rowind, const std::vector<int>& col, int nrhs);

    // Destructor
    virtual ~LapackQRDenseInternal();
    
    // Initialize the solver
    virtual void init();

    // Prepare the solution of the linear system
    virtual void prepare();
    
    // Solve the system of equations
    virtual void solve();

  protected:

    // Matrix
    std::vector<double> mat_;
    
    // The scalar factors of the elementary reflectors
    std::vector<double> tau_; 
    
    // qr work array
    std::vector<double> work_; 
    
};



} // namespace CasADi

#endif //LAPACK_QR_DENSE_HPP

