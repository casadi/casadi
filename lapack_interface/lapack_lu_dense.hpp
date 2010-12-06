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

#ifndef LAPACK_LU_DENSE_HPP
#define LAPACK_LU_DENSE_HPP

#include "casadi/fx/linear_solver.hpp"

namespace CasADi{
  
/** \brief  Forward declaration of internal class */
class LapackLUDenseInternal;

/** \brief  Public class */
class LapackLUDense : public LinearSolver{
public:

  /// Default (empty) constructor
  LapackLUDense();
  
  /// Create a linear solver given a sparsity pattern
  LapackLUDense(int nrow, int ncol, const std::vector<int>& rowind, const std::vector<int>& col, int nrhs=1);
    
  /// Access functions of the node
  LapackLUDenseInternal* operator->();
  const LapackLUDenseInternal* operator->() const;
};

/// LU-Factorize dense matrix (lapack)
extern "C" void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

/// Solve a system of equation using an LU-factorized matrix (lapack)
extern "C" void dgetrs_(char* trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

/// Calculate row and column scalings
extern "C" void dgeequ_(int *m, int *n, double *a, int *lda, double *r, double *c, double *rowcnd, double *colcnd, double *amax, int *info);

/// Equilibriate the system
extern "C" void dlaqge_(int *m, int *n, double *a, int *lda, double *r, double *c, double *rowcnd, double *colcnd, double *amax, char *equed );

/// Internal class
class LapackLUDenseInternal : public LinearSolverInternal{
  public:
    // Create a linear solver given a sparsity pattern and a number of right hand sides
    LapackLUDenseInternal(int nrow, int ncol, const std::vector<int>& rowind, const std::vector<int>& col, int nrhs);

    // Destructor
    virtual ~LapackLUDenseInternal();
    
    // Initialize the solver
    virtual void init();

    // Prepare the solution of the linear system
    virtual void prepare();
    
    // Solve the system of equations
    virtual void solve();

  protected:

    // Matrix
    std::vector<double> mat_;
    
    // Pivoting elements
    std::vector<int> ipiv_;

    // Row and column scaling
    std::vector<double> r_, c_;
    
    // Type of scaling during the last equilibration
    char equed_;
    
    // Equilibriate?
    bool equilibriate_;

    // Allow the equilibration to fail
    bool allow_equilibration_failure_;
    
    // Factorization sucessful?
    bool factorization_sucessful_;
    
};



} // namespace CasADi

#endif //LAPACK_LU_DENSE_HPP

