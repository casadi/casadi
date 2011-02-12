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

#ifndef OCP_SOLVER_INTERNAL_HPP
#define OCP_SOLVER_INTERNAL_HPP

#include <vector>
#include "ocp_solver.hpp"
#include "fx_internal.hpp"

namespace CasADi{
 
  /** \brief  Internal node class for OCPSolver
  \author Joel Andersson 
  \date 2010
*/
class OCPSolverInternal : public FXInternal{
  friend class OCPSolver;
  public:
  
    /// Constructor
    explicit OCPSolverInternal(const FX& ffcn, const FX& mfcn, const FX& cfcn, const FX& rfcn);

    /// Destructor
    virtual ~OCPSolverInternal();
    
    /// Initialize
    virtual void init();

    // Discrete time dynamics
    FX ffcn_;
    
    // Mayer term
    FX mfcn_;
    
    // Path constraints
    FX cfcn_;
    
    // Initial value constraints
    FX rfcn_;

    // Number of grid points
    int nk_;

    // Number of differential states
    int nx_;

    // Number of algebraic states
    int nz_;

    // Number of parameters
    int np_;
    
    // Number of controls
    int nu_;
    
    // Number of point constraints
    int nh_;
    
    // Number of point coupling constraints
    int ng_;
    
    /// Cost functions
/*    std::vector<FX> L_;
    
    /// Dynamic constraint
    std::vector<FX> F_;
    
    /// Point constraints
    std::vector<FX> H_;
    
    /// Coupling constraints
    std::vector<FX> G_;*/
};



} // namespace CasADi


#endif // OCP_SOLVER_INTERNAL_HPP

