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

#ifndef SDP_SOLVER_INTERNAL_HPP
#define SDP_SOLVER_INTERNAL_HPP

#include "sdp_solver.hpp"
#include "fx_internal.hpp"

namespace CasADi{

/// Internal class
class SDPSolverInternal : public FXInternal{
  public:
    // Constructor
    SDPSolverInternal();
        
    // Constructor
    SDPSolverInternal( const CRSSparsity &C, const CRSSparsity &A);
    
    // Destructor
    virtual ~SDPSolverInternal() = 0;
    
    // Initialize
    virtual void init();
    
    // Solve the system of equations
    virtual void evaluate(int nfdir, int nadir);
    
    // Solve the system of equations
    virtual void solve();
    
  protected:
    
    /// Size of decision variable matrix
    int n_;
    
    /// The number of matrix constraints
    int m_;
    
    /// A mapping from A to Ai
    FX mapping_;
    
    /** \brief Indicate if dual should be allocated and calculated.
    * You may want to avoid calculating the dual variable for problems with n large, as the dual (n x n) is always dense.
    */
    bool calc_dual_;
    
    /* \brief indicate if P part of primal solution should be allocated and calculated
    * You may want to avoid calculating this variable for problems with n large, as is always dense (n x n).
    */
    bool calc_p_;
};


} // namespace CasADi

#endif //SDP_SOLVER_INTERNAL_HPP

