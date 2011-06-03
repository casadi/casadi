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

#include "fx.hpp"
#include "linear_solver.hpp"

namespace CasADi{

/// Input arguments of an ODE/DAE function
enum DAEInput{
  /** Time. (1-by-1) */
  DAE_T,
  /** State vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES */
  DAE_Y,
  /** Parameter vector (matrix). */
  DAE_P,
  /** State derivative vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES */
  DAE_YDOT,
  /** Number of arguments. */
  DAE_NUM_IN
};

/// Output arguments of an ODE/DAE residual function
enum DAEOutput{
  /** Right hand side of ODE. Should have same amount of non-zeros as ODEInput:ODE_Y */
  DAE_RES,
  /** Number of arguments. */
  DAE_NUM_OUT
};

/// Input arguments of an integrator
enum IntegratorInput{
  /** Differential or algebraic state at t0  (dimension nx-by-1) */
  INTEGRATOR_X0, 
  /** Parameters p  (dimension np-by-1) */
  INTEGRATOR_P,  
  /** State derivative at t0  (dimension nx-by-1)
  * This input may be changed during an IDASIntegrator::evaluate()
  */
  INTEGRATOR_XP0, 
  /** Number of input arguments of an integrator */
  INTEGRATOR_NUM_IN};

/// Output arguments of an integrator
enum IntegratorOutput{
 /**  State at tf */
 INTEGRATOR_XF, 
 /**  State derivative at tf */
 INTEGRATOR_XPF, 
  /** Number of output arguments of an integrator */
 INTEGRATOR_NUM_OUT
};

/// Forward declaration of internal class
class IntegratorInternal;

/** Integrator abstract base class
  An "integrator" is a function that solves an initial value problem (IVP) of the generic form:
  
  F(t,x,der(x),z,p) == 0
  x(t0) = x0
  over a time interval [t0, tf].

  It has (currently) 6 inputs, initial time, final time, initial state (vector-valued), parameter (vector-valued), as well as guesses for the initial 
  state derivative and algebraic variables. The latter two are only relevant for implicit integrators.
  
  In addition to this, the integrator provides some additional functionality, such as getting the value of the state and/or sensitivities at certain time points.
  Controls are assumed to be parametrized at this point. 
    
  The class does not specify how the function F above should be represented, nor the method used for the integration, but assumes that it steps forward in time 
  (ruling out collocation in particular). The actual form of the ODE/DAE is defined in the derived classes.
    
  inputs:
  0: State at t0  (dimension nx-by-1)
  1: Parameter  (dimension np-by-1)
  2: State derivative at t0  (dimension nx-by-1)
  
  outputs:
  0: State at tf
  1: State derivative at tf


  \author Joel Andersson
  \date 2010
*/

// grep "addOption" integrator_internal.cpp | perl -pe 's/addOption\((.*?),(.*?),(.*?)\);(.*\/\/ (.*))?/* \1 \2 \3 ...  \5\\n/'

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
  void reset(int fsens_order=0, int asens_order=0);

  /// Integrate until a specified time point 
  void integrate(double t_out);

  /// Reset the solver of the adjoint problem and take time to tf
  void resetAdj();

  /// Integrate backwards in time until a specified time point
  void integrateAdj(double t_out);

  /// Set initial time
  void setInitialTime(double t0);

  /// Set final time
  void setFinalTime(double tf);

  /// Set a stop time for the forward integration
  void setStopTime(double tf);

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;

  /// Set linear solver
  void setLinearSolver(const LinearSolver& linsol, const FX& jac=FX());
  
  /// Get the Jacobian
  FX getJacobian();
  
  /// Get the Linear solver
  LinearSolver getLinearSolver();

  /// Generate a new integrator integrating the forward sensitivity augmented ODE/DAE
  Integrator jac(int iind=0, int oind=0);
};


} // namespace CasADi

#endif //INTEGRATOR_HPP
