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

#ifndef OCP_SOLVER_HPP
#define OCP_SOLVER_HPP

#include <vector>

#include "fx.hpp"

namespace CasADi{

  /// Input arguments of an OCP Solver \n
  ///   ns: Number of shooting nodes: from option number_of_grid_points\n
  ///   nx: Number of differential states: from ffcn.input(INTEGRATOR_X0).size() \n
  ///   nu: Number of controls: from ffcn.input(INTEGRATOR_P).size() - np \n
  ///   np: Number of parameters: from option number_of_parameters\n
  ///   nh: Number of point constraints: from cfcn.input(0).size()
  enum OCPInput{
    /// Time grid: ((ns+1) x 1) - default: linspace(0,t_final,ns+1) 
    OCP_T,   
    /// States lower bounds (nx x (ns+1))
    OCP_LBX,
    /// States upper bounds (nx x (ns+1))
    OCP_UBX,
    /// States initial guess (nx x (ns+1))
    OCP_X_INIT,
    /// States deriatives lower bounds (nx x (ns+1))
    OCP_LBXP,
    /// States deriatives upper bounds (nx x (ns+1))
    OCP_UBXP,
    /// States deriatives initial guess (nx x (ns+1))
    OCP_XP_INIT,
    /// Controls lower bounds (nu x ns)
    OCP_LBU,
    /// Controls upper bounds (nu x ns)
    OCP_UBU,
    /// Controls initial guess (nu x ns)
    OCP_U_INIT,     
    /// Parameters lower bounds (np x 1)
    OCP_LBP,
    /// Parameters upper bounds (np x 1)
    OCP_UBP,
    /// Parameters initial guess (np x 1)
    OCP_P_INIT,
    /// Point constraint lower bound (nh x (ns+1))
    OCP_LBH,
    /// Point constraint upper bound (nh x (ns+1))
    OCP_UBH,
    /// Bounds for the coupling constraints
    OCP_LBG, OCP_UBG,
    /// Number of inputs to an OCP solver
    OCP_NUM_IN
  };
  
  /// Outputs arguments of an OCP Solver
  enum OCPOutput{
    /// Optimal state trajectory
    OCP_X_OPT, 
    /// Optimal control trajectory
    OCP_U_OPT, 
    /// Optimal state derivative trajectory
    OCP_XP_OPT, 
    /// Number of outputs to an OCP solver
    OCP_NUM_OUT
  };
    
  // Forward declaration of internal class
  class OCPSolverInternal;

  /** \brief Base class for OCP solvers
   *
   *
   * OCPSolver is an CasADi::FX mapping from CasADi::OCPInput to CasADi::OCPOutput
   *
      \author Joel Andersson
      \date 2011
  */ 
  class OCPSolver : public FX{
  public:

    /// Default constructor
    OCPSolver();

    /// Access functions of the node
    OCPSolverInternal* operator->();

    /// Const access functions of the node
    const OCPSolverInternal* operator->() const;
    
  };

} // namespace CasADi


#endif // OCP_SOLVER_HPP

