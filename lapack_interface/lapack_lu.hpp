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

#ifndef LAPACK_LU_HPP
#define LAPACK_LU_HPP

#include "casadi/fx/linear_solver.hpp"

namespace CasADi{
  
/** \brief  Forward declaration of internal class */
class LapackLUInternal;

/** \brief  Public class */
class LapackLU : public LinearSolver{
public:

  /// Default (empty) constructor
  LapackLU();
  
  /// Create a linear solver given a sparsity pattern
  LapackLU(int nrow, int ncol, int nrhs=1);
    
  /// Access functions of the node
  LapackLUInternal* operator->();
  const LapackLUInternal* operator->() const;
};

/// LU-Factorize dense matrix (lapack)
extern "C" void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

/// Solve a system of equation using an LU-factorized matrix (lapack)
extern "C" void dgetrs_(char* trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

/// Internal class
class LapackLUInternal : public LinearSolverInternal{
  public:
    // Create a linear solver given a sparsity pattern and a number of right hand sides
    LapackLUInternal(int nrow, int ncol, int nrhs);

    // Destructor
    virtual ~LapackLUInternal();
    
    // Initialize the solver
    virtual void init();

    // Prepare the solution of the linear system
    virtual void prepare();
    
    // Solve the system of equations
    virtual void solve();

  protected:
    int nrow_, ncol_;
    int nrhs_;
};



} // namespace CasADi

#endif //LAPACK_LU_HPP

