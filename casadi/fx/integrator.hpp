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

#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include "../expression_tools.hpp"
#include "fx.hpp"
#include "linear_solver.hpp"
#include "integrator_jacobian.hpp"

namespace CasADi{

/// Input arguments of an integrator
enum IntegratorInput{INTEGRATOR_T0, INTEGRATOR_TF, INTEGRATOR_X0, INTEGRATOR_P, INTEGRATOR_XP0, INTEGRATOR_NUM_IN};

// Input arguments of an integrator (new)
//enum IntegratorInput{INTEGRATOR_X0, INTEGRATOR_P, INTEGRATOR_NUM_IN};

/// Output arguments of an integrator
enum IntegratorOutput{INTEGRATOR_XF, INTEGRATOR_XPF, INTEGRATOR_NUM_OUT};

// Output arguments of an integrator (new)
//enum IntegratorOutput{INTEGRATOR_XF, INTEGRATOR_NUM_OUT};

/// Forward declaration of internal class
class IntegratorInternal;

/** Integrator abstract base class
  An "integrator" is a function that solves an initial value problem (IVP) of the generic form:
  
  F(x,der(x),p,t) == 0
  x(t0) = x0

  It has 4 inputs, initial time, final time, initial state (vector-valued) and parameter (vector-valued) and one output, the state at the final time. 
  In addition to this, the integrator provides some additional functionality, such as getting the value of the state and/or sensitivities at certain time points.
  Controls are assumed to be parametrized at this point. 
    
  The class does not specify how the function F above should be represented, nor the method used for the integration, but assumes that it steps forward in time (ruling out collocation in particular). 
  The actual form of the ODE/DAE is defined in the derived classes.
    
  inputs:
  0: Initial time t0 (dimension 1-by-1)
  1: Final time tf (dimension 1-by-1)
  2: State at t0  (dimension nx-by-1)
  3: p  (dimension np-by-1)
  
  outputs:
  0: y(tf)
  
  \author Joel Andersson
  \date 2010
*/
class Integrator : public FX{
public:
  /// Default constructor
  Integrator();

  /// Clone
  Integrator clone() const;
  
  /// Print solver statistics
  void printStats(std::ostream &stream=std::cout) const;
  
  /// Access functions of the node
  IntegratorInternal* operator->();

  /// Access functions of the node
  const IntegratorInternal* operator->() const;
  
  /// Reset the solver and bring the time back to t0 
  void reset(int fsens_order, int asens_order=0);

  /// Integrate until a specified time point 
  void integrate(double t_out);

  /// Reset the solver of the adjoint problem and take time to tf
  void resetAdj();

  /// Integrate backwards in time until a specified time point
  void integrateAdj(double t_out);

  /// Set a stop time for the forward integration
  void setStopTime(double tf);

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;

  /// Set linear solver
  void setLinearSolver(const LinearSolver& linsol, const FX& jac=FX());
  
  /// Jacobian of output oind with respect to input iind
  IntegratorJacobian jacobian(int iind=0, int oind=0);

  /// Generate a new integrator integrating the forward sensitivity augmented ODE/DAE
  Integrator jac(int iind=0, int oind=0);
};


} // namespace CasADi

#endif //INTEGRATOR_HPP
