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

#include "function.hpp"

namespace CasADi{


/// Input arguments of a Mayer Term \n [mayerIn]
///   nx: Number of states: from ffcn.input(INTEGRATOR_X0).size() \n
///   np: Number of parameters: from option number_of_parameters\n
enum MayerInput{
  /// States at the end of integration (nx x 1) [x]
  MAYER_X,
  /// Problem parameters (np x 1) [p]
  MAYER_P,
  /// Number of inputs to a mayer term
  MAYER_NUM_IN
};

/// Input arguments of an OCP Solver \n [ocpIn]
///   ns: Number of shooting nodes: from option number_of_grid_points\n
///   nx: Number of states: from ffcn.input(INTEGRATOR_X0).size() \n
///   nc: Number of constants duting intergation: ffcn.input(INTEGRATOR_P).size()
///   nu: Number of controls: from nc - np \n
///   np: Number of parameters: from option number_of_parameters\n
///   nh: Number of point constraints: from cfcn.input(0).size()
enum OCPInput{
  /// States lower bounds (nx x (ns+1)) [lbx]
  OCP_LBX,
  /// States upper bounds (nx x (ns+1)) [ubx]
  OCP_UBX,
  /// States initial guess (nx x (ns+1)) [x_init]
  OCP_X_INIT,
  /// Controls lower bounds (nu x ns) [lbu]
  OCP_LBU,
  /// Controls upper bounds (nu x ns) [ubu]
  OCP_UBU,
  /// Controls initial guess (nu x ns) [u_init]
  OCP_U_INIT,     
  /// Parameters lower bounds (np x 1) [lbp]
  OCP_LBP,
  /// Parameters upper bounds (np x 1) [ubp]
  OCP_UBP,
  /// Parameters initial guess (np x 1) [p_init]
  OCP_P_INIT,
  /// Point constraint lower bound (nh x (ns+1)) [lbh]
  OCP_LBH,
  /// Point constraint upper bound (nh x (ns+1)) [ubh]
  OCP_UBH,
  /// Lower bound for the coupling constraints [lbg]
  OCP_LBG,
  /// Upper bound for the coupling constraints [ubg]
  OCP_UBG,
  /// Number of inputs to an OCP solver
  OCP_NUM_IN
};

/// Output arguments of an OCP Solver [ocpOut]
enum OCPOutput{
  /// Optimal state trajectory [x_opt]
  OCP_X_OPT, 
  /// Optimal control trajectory [u_opt]
  OCP_U_OPT, 
  /// Optimal parameters [p_opt]
  OCP_P_OPT, 
  /// Objective/cost function for optimal solution (1 x 1) [cost]
  OCP_COST,
  /// Number of outputs to an OCP solver
  OCP_NUM_OUT
};
    
  // Forward declaration of internal class
  class OCPSolverInternal;

  /** \brief Base class for OCP solvers
   *
   *
      \author Joel Andersson
      \date 2011-2013
  */ 
  class OCPSolver : public Function{
  public:

    /// Default constructor
    OCPSolver();

    /// Access functions of the node
    OCPSolverInternal* operator->();

    /// Const access functions of the node
    const OCPSolverInternal* operator->() const;
    
    // Access the underlying ffcn
    Function getFfcn() const;
    
    // Access the underlying mfcn
    Function getMfcn() const;
    
    // Access the underlying cfcn
    Function getCfcn() const;
    
    // Access the underlying rfcn
    Function getRfcn() const;
    
    
  };

} // namespace CasADi


#endif // OCP_SOLVER_HPP

