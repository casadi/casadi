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
#include "integrator_jacobian.hpp"

namespace CasADi{

/// Input arguments of an integrator
enum IntegratorInput{
  /** Initial time t0 (dimension 1-by-1) */
  INTEGRATOR_T0, 
  /** Final time tf (dimension 1-by-1) */
  INTEGRATOR_TF,
  /** Differential State at t0  (dimension nx-by-1) */
  INTEGRATOR_X0, 
  /** Parameters p  (dimension np-by-1) */
  INTEGRATOR_P,  
  /** Differential state derivative at t0  (dimension nx-by-1) */
  INTEGRATOR_XP0, 
  /** Algebraic state at t0  (dimension nz-by-1) */
  INTEGRATOR_Z0, 
  INTEGRATOR_NUM_IN};

/// Output arguments of an integrator
enum IntegratorOutput{
 /**  Differential state at tf */
 INTEGRATOR_XF, 
 /**  Differential state derivative at tf */
 INTEGRATOR_XPF, 
  /** Algebraic state at tf*/
 INTEGRATOR_ZF, 
 INTEGRATOR_NUM_OUT};

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
  0: Initial time t0 (dimension 1-by-1)
  1: Final time tf (dimension 1-by-1)
  2: Differential state at t0  (dimension nx-by-1)
  3: p  (dimension np-by-1)
  4: Differential state derivative at t0  (dimension nx-by-1)
  5: Algebraic state at t0  (dimension nz-by-1)
  
  outputs:
  0: Differential state at tf
  1: Differential state derivative at tf
  2: Algebraic state at tf
  
  Options:
  
  * "max_num_steps"                OT_INTEGER  10000 ...  maximum number of steps\n
  * "reltol"                       OT_REAL     1e-6 ...  relative tolerence for the IVP solution\n
  * "abstol"                       OT_REAL     1e-8 ...  absolute tolerence  for the IVP solution\n
  * "upper_bandwidth"              OT_INTEGER   Option() ...  upper band-width of banded jacobians\n
  * "lower_bandwidth"              OT_INTEGER   Option() ...  lower band-width of banded jacobians\n
  * "linear_solver"                OT_STRING  "dense" ...  "dense", "banded" or "iterative"\n
  * "iterative_solver"             OT_STRING  "gmres" ...  "gmres", "bcgstab", "tfqmr"\n
  * "pretype"                      OT_STRING  "none" ...  "none", "left", "right", "both"\n
  * "exact_jacobian"               OT_BOOLEAN   false ...  \n
  * "max_krylov"                   OT_INTEGER   10 ...  maximum krylov subspace size\n
  * "is_differential"              OT_INTEGERVECTOR   Option() ...  \n
  * "sensitivity_method"           OT_STRING   "simultaneous" ...  "simultaneous" or "staggered"\n
  * "max_multistep_order"          OT_INTEGER  5 ...  \n
  * "use_preconditioner"           OT_BOOLEAN  false ...  precondition an iterative solver\n
  * "stop_at_end"                  OT_BOOLEAN  false ...  Stop the integrator at the end of the interval\n
  * "jacmap"                       OT_INTEGERVECTOR  Option() ...  if the integrator is the Jacobian of another integrator, this option will contain the mapping between the states\n
  * "jacinit"                      OT_REALVECTOR  Option() ...  initial values to the forward sensitivities\n
  * "nrhs"                         OT_INTEGER  1 ...  number of right hand sides\n
  * "quad_err_con"                 OT_BOOLEAN false ...  should the quadratures affect the step size control\n
  * "fsens_err_con"                OT_INTEGER  false ...  include the forward sensitivities in all error controls\n
  * "finite_difference_fsens"      OT_BOOLEAN  false ...  use finite differences to approximate the forward sensitivity equations (if AD is not available)\n
  * "fsens_reltol"                 OT_REAL     Option() ...  relative tolerence for the forward sensitivity solution [default: equal to reltol]\n
  * "fsens_abstol"                 OT_REAL     Option() ...  absolute tolerence for the forward sensitivity solution [default: equal to abstol]\n
  * "fsens_scaling_factors"        OT_REALVECTOR  Option() ...  scaling factor for the components if finite differences is used\n
  * "fsens_sensitiviy_parameters"  OT_INTEGERVECTOR  Option() ...  specifies which components will be used when estimating the sensitivity equations\n
  * "steps_per_checkpoint"         OT_INTEGER 20 ...  number of steps between two consecutive checkpoints\n
  * "interpolation_type"           OT_STRING "hermite" ...  type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")\n
  * "asens_upper_bandwidth"        OT_INTEGER   Option() ...  upper band-width of banded jacobians\n
  * "asens_lower_bandwidth"        OT_INTEGER   Option() ...  lower band-width of banded jacobians\n
  * "asens_linear_solver"          OT_STRING  "dense" ...  "dense", "banded" or "iterative"\n
  * "asens_iterative_solver"       OT_STRING  "gmres" ...  "gmres", "bcgstab", "tfqmr"\n
  * "asens_pretype"                OT_STRING  "none" ...  "none", "left", "right", "both"\n
  * "asens_max_krylov"             OT_INTEGER   10 ...  maximum krylov subspace size\n
  * "asens_reltol"                      OT_REAL     Option() ...  relative tolerence for the adjoint sensitivity solution [default: equal to reltol]\n
  * "asens_abstol"                      OT_REAL     Option() ...  absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]\n

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
  
  /** \brief Jacobian of output oind with respect to input iind
  *
  *  Jacobian is implemented by augmenting the system with equations for the time evolution of sensitivity.
  *  The result is a mapping from CasADi::IntegratorInput to CasADi::IntegratorJacobianOutput.
  *
  */
  IntegratorJacobian jacobian(int iind=0, int oind=0);

  /// Generate a new integrator integrating the forward sensitivity augmented ODE/DAE
  Integrator jac(int iind=0, int oind=0);
};


} // namespace CasADi

#endif //INTEGRATOR_HPP
