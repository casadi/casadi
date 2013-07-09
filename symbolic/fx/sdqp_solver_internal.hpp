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

#ifndef SDQP_SOLVER_INTERNAL_HPP
#define SDQP_SOLVER_INTERNAL_HPP

#include "sdqp_solver.hpp"
#include "fx_internal.hpp"

namespace CasADi{

/// Internal class
class SDQPSolverInternal : public FXInternal{
  public:

    // Constructor
    SDQPSolverInternal(const std::vector<CRSSparsity>& st);
    
    // Destructor
    virtual ~SDQPSolverInternal() = 0;
    
    // Initialize
    virtual void init();
    
    // Solve the system of equations
    virtual void evaluate(int nfdir, int nadir);
    
    // Solve the system of equations
    virtual void solve();
    
    /// Set options that make the SDQP solver more suitable for solving SOCPs
    virtual void setSOCQPOptions() { };
    
    /// Print out problem statement for debugging
    void printProblem(std::ostream &stream=std::cout) const;
    
  protected:
  
    /// Problem structure
    std::vector<CRSSparsity> st_;
    
    /// Size of decision variable vector
    int n_;
    
    /// The number of matrices F
    int m_;
    
    /// Number of identified blocks
    int nb_;
    
    /// Number of linear constraints
    int nc_;

};


} // namespace CasADi

#endif //SDQP_SOLVER_INTERNAL_HPP

