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

#ifndef QP_SOLVER_INTERNAL_HPP
#define QP_SOLVER_INTERNAL_HPP

#include "qp_solver.hpp"
#include "fx_internal.hpp"

namespace CasADi{

/// Internal class
class QPSolverInternal : public FXInternal{
  public:

    // Constructor
    QPSolverInternal(const std::vector<CRSSparsity> &st);
    
    // Destructor
    virtual ~QPSolverInternal() = 0;
    
    // Initialize
    virtual void init();
    
    // Solve the system of equations
    virtual void evaluate(int nfdir, int nadir);
    
    // Solve the system of equations
    virtual void solve();
    
    /// Set options that make the QP solver more suitable for solving LPs
    virtual void setLPOptions() { };
    
  protected:

    /// Problem structure
    std::vector<CRSSparsity> st_;
    
    /// Number of decision variables
    int nx_;
    
    /// The number of constraints (counting both equality and inequality) == A.size1()
    int nc_; 
};


} // namespace CasADi

#endif //QP_SOLVER_INTERNAL_HPP

