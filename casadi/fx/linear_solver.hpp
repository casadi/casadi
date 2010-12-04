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

#ifndef LINEAR_SOLVER_HPP
#define LINEAR_SOLVER_HPP

#include "fx.hpp"

namespace CasADi{
  
/// Forward declaration of internal class
class LinearSolverInternal;

/// Public class
class LinearSolver : public FX{
public:
  
  /// Access functions of the node
  LinearSolverInternal* operator->();
  const LinearSolverInternal* operator->() const;
  
  /// Factorize the matrix
  void prepare();
    
  /// Solve the system of equations
  void solve();
};

/// Internal class
class LinearSolverInternal : public FXNode{
  public:
    // Constructor
    LinearSolverInternal(int nrow, int ncol, const std::vector<int>& rowind, const std::vector<int>& col, int nrhs);
    
    // Destructor
    virtual ~LinearSolverInternal() = 0;
    
    // Solve the system of equations
    virtual void evaluate(int fsens_order, int asens_order);

    // Prepare the factorization
    virtual void prepare() = 0;
    
    // Solve the system of equations
    virtual void solve() = 0;
    
    // Sparsity in CRS format
    int nrow_, ncol_;
    std::vector<int> rowind_, col_;
    int nrhs_;
};


} // namespace CasADi

#endif //LINEAR_SOLVER_HPP

