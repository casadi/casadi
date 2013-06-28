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

#ifndef LINEAR_SOLVER_INTERNAL_HPP
#define LINEAR_SOLVER_INTERNAL_HPP

#include "linear_solver.hpp"
#include "fx_internal.hpp"

namespace CasADi{
  

  /** Internal class
      @copydoc LinearSolver_doc
  */
  class LinearSolverInternal : public FXInternal{
  public:
    // Constructor
    LinearSolverInternal(const CRSSparsity& sparsity, int nrhs);
        
    // Destructor
    virtual ~LinearSolverInternal() = 0;
    
    // Initialize
    virtual void init();
    
    // Solve the system of equations
    virtual void evaluate(int nfdir, int nadir);

    // Prepare the factorization
    virtual void prepare() = 0;

    // Solve the system of equations, using internal vector
    virtual void solve();

    // Solve the system of equations
    virtual void solve(double* x, int nrhs, bool transpose) = 0;

    // Is prepared
    bool prepared_;

    // Get sparsity pattern
    int nrow() const{ return input(LINSOL_A).size1();}
    int ncol() const{ return input(LINSOL_A).size2();}
    int nnz() const{ return input(LINSOL_A).size();}
    const std::vector<int>& col() const{ return input(LINSOL_A).col();}
    const std::vector<int>& rowind() const{ return input(LINSOL_A).rowind();}    
  };


} // namespace CasADi

#endif //LINEAR_SOLVER_INTERNAL_HPP

