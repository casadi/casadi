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

#include "integrator_internal.hpp"
#include <cassert>
#include "../stl_vector_tools.hpp"
#include "jacobian.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "mx_function.hpp"

using namespace std;
namespace CasADi{

IntegratorInternal::IntegratorInternal(){
  // set default options
  setOption("name","unnamed_integrator"); // name of the function
  
  // IVP solution
  addOption("max_num_steps",               OT_INTEGER, 10000); // maximum number of steps
  addOption("reltol",                      OT_REAL,    1e-6); // relative tolerence for the IVP solution
  addOption("abstol",                      OT_REAL,    1e-8); // absolute tolerence  for the IVP solution
  addOption("upper_bandwidth",             OT_INTEGER); // upper band-width of banded jacobians
  addOption("lower_bandwidth",             OT_INTEGER); // lower band-width of banded jacobians
  addOption("linear_solver",               OT_STRING, "dense"); // "dense", "banded" or "iterative"
  addOption("iterative_solver",            OT_STRING, "gmres"); // "gmres", "bcgstab", "tfqmr"
  addOption("pretype",                     OT_STRING, "none"); // "none", "left", "right", "both"
  addOption("exact_jacobian",              OT_BOOLEAN,  false);
  addOption("max_krylov",                  OT_INTEGER,  10);        // maximum krylov subspace size
  addOption("is_differential",             OT_INTEGERVECTOR);
  addOption("sensitivity_method",          OT_STRING,  "simultaneous"); // "simultaneous" or "staggered"
  addOption("max_multistep_order",         OT_INTEGER, 5);
  addOption("use_preconditioner",          OT_BOOLEAN, false); // precondition an iterative solver
  addOption("stop_at_end",                 OT_BOOLEAN, false); // Stop the integrator at the end of the interval
  addOption("nrhs",                        OT_INTEGER, 1); // number of right hand sides
  addOption("t0",                          OT_REAL, 0.0); // start of the integration
  addOption("tf",                          OT_REAL, 1.0); // end of the integration

  // Quadratures
  addOption("quad_err_con",                OT_BOOLEAN,false); // should the quadratures affect the step size control

  // Forward sensitivity problem
  addOption("fsens_err_con",               OT_INTEGER, false); // include the forward sensitivities in all error controls
  addOption("finite_difference_fsens",     OT_BOOLEAN, false); // use finite differences to approximate the forward sensitivity equations (if AD is not available)
  addOption("fsens_reltol",                OT_REAL); // relative tolerence for the forward sensitivity solution [default: equal to reltol]
  addOption("fsens_abstol",                OT_REAL); // absolute tolerence for the forward sensitivity solution [default: equal to abstol]
  addOption("fsens_scaling_factors",       OT_REALVECTOR); // scaling factor for the components if finite differences is used
  addOption("fsens_sensitiviy_parameters", OT_INTEGERVECTOR); // specifies which components will be used when estimating the sensitivity equations

  // Adjoint sensivity problem
  addOption("steps_per_checkpoint",        OT_INTEGER,20); // number of steps between two consecutive checkpoints
  addOption("interpolation_type",          OT_STRING,"hermite"); // type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")
  addOption("asens_upper_bandwidth",       OT_INTEGER); // upper band-width of banded jacobians
  addOption("asens_lower_bandwidth",       OT_INTEGER); // lower band-width of banded jacobians
  addOption("asens_linear_solver",         OT_STRING, "dense"); // "dense", "banded" or "iterative"
  addOption("asens_iterative_solver",      OT_STRING, "gmres"); // "gmres", "bcgstab", "tfqmr"
  addOption("asens_pretype",               OT_STRING, "none"); // "none", "left", "right", "both"
  addOption("asens_max_krylov",            OT_INTEGER,  10);        // maximum krylov subspace size
  addOption("asens_reltol",                OT_REAL); // relative tolerence for the adjoint sensitivity solution [default: equal to reltol]
  addOption("asens_abstol",                OT_REAL); // absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]
  
  nx_ = 0;
  np_ = 0;
}

IntegratorInternal::~IntegratorInternal(){ 
}

void IntegratorInternal::setDimensions(int nx, int np){
  nx_ = nx;
  np_ = np;
  
  // Allocate space for inputs
  input_.resize(INTEGRATOR_NUM_IN);
  input(INTEGRATOR_X0)  = DMatrix(nx_,1,0); // initial state value
  input(INTEGRATOR_XP0) = DMatrix(nx_,1,0); // initial state derivative value
  input(INTEGRATOR_P)   = DMatrix(np_,1,0); // parameter
  
  // Allocate space for outputs
  output_.resize(INTEGRATOR_NUM_OUT);
  output(INTEGRATOR_XF) = DMatrix(nx_,1,0);
  output(INTEGRATOR_XPF)= DMatrix(nx_,1,0);
}

void IntegratorInternal::evaluate(int nfdir, int nadir){
  
  // Reset solver
  reset(nfdir>0, nadir>0);

  // Set the stop time of the integration -- don't integrate past this point
  if(stop_at_end_) setStopTime(tf_);

  // Advance solution in time
  integrate(tf_);

  if(nadir>0){
    // Re-initialize backward problem
    resetAdj();

    // Integrate backwards to the beginning
    integrateAdj(t0_);
  }
}

void IntegratorInternal::init(){
  // Call the base class method
  FXInternal::init();

  // read options
  exact_jacobian_ = getOption("exact_jacobian");
  abstol_ = getOption("abstol"); // TODO: change to vector tolerences
  reltol_ = getOption("reltol");
  max_num_steps_ = getOption("max_num_steps");
  finite_difference_fsens_ = getOption("finite_difference_fsens").toInt();
  fsens_abstol_ = hasSetOption("fsens_abstol") ? double(getOption("fsens_abstol")) : abstol_;
  fsens_reltol_ = hasSetOption("fsens_reltol") ? double(getOption("fsens_reltol")) : reltol_;
  asens_abstol_ = hasSetOption("asens_abstol") ? double(getOption("asens_abstol")) : abstol_;
  asens_reltol_ = hasSetOption("asens_reltol") ? double(getOption("asens_reltol")) : reltol_;
  stop_at_end_ = getOption("stop_at_end");
  nrhs_ = getOption("nrhs");
  
  // Give an intial value for the time horizon
  setInitialTime(getOption("t0"));
  setFinalTime(getOption("tf"));
}

FX IntegratorInternal::jacobian(int iind, int oind){
  if(symbjac()){
  
  // Create a new integrator for the forward sensitivity equations
  Integrator fwdint = jac(iind,oind);
  fwdint.init();

  // Number of sensitivities
  int ns;
  switch(iind){
    case INTEGRATOR_P: ns = np_; break;
    case INTEGRATOR_X0: ns = nx_; break;
    default: casadi_assert_message(false,"iind must be INTEGRATOR_P or INTEGRATOR_X0");
  }
  
  // Symbolic input of the Jacobian
  vector<MX> jac_in = symbolicInput();
  
  // Input to the augmented integrator
  vector<MX> fwdint_in(INTEGRATOR_NUM_IN);
  
  // Pass parameters without change
  fwdint_in[INTEGRATOR_P] = jac_in[INTEGRATOR_P];
  
  // Get the state
  MX x0 = jac_in[INTEGRATOR_X0];
  MX xp0 = jac_in[INTEGRATOR_XP0];
  
  // Separate the quadrature states from the rest of the states
  MX y0 = x0[range(ny_)];
  MX q0 = x0[range(ny_,nx_)];
  MX yp0 = xp0[range(ny_)];
  MX qp0 = xp0[range(ny_,nx_)];
  
  // Initial condition for the sensitivitiy equations
  DMatrix y0_sens(ns*ny_,1,0);
  DMatrix q0_sens(ns*nq_,1,0);
  
  if(iind == INTEGRATOR_X0){
    for(int i=0; i<ny_; ++i){
      y0_sens[i + i*ns] = 1;
    }
    
    for(int i=0; i<nq_; ++i){
      q0_sens[ny_ + i + i*ns] = 1;
    }
  }
  
  // Augmented initial condition
  MX y0_aug = vertcat(y0,MX(y0_sens));
  MX q0_aug = vertcat(q0,MX(q0_sens));
  MX yp0_aug = vertcat(yp0,MX::zeros(y0_sens.sparsity()));
  MX qp0_aug = vertcat(qp0,MX::zeros(q0_sens.sparsity()));

  // Finally, we are ready to pass the initial condition for the state and state derivative
  fwdint_in[INTEGRATOR_X0] = vertcat(y0_aug,q0_aug);
  fwdint_in[INTEGRATOR_XP0] = vertcat(yp0_aug,qp0_aug);
  
  // Call the integrator with the constructed input (in fact, create a call node)
  vector<MX> fwdint_out = fwdint.call(fwdint_in);
  MX xf_aug = fwdint_out[INTEGRATOR_XF];
  MX xpf_aug = fwdint_out[INTEGRATOR_XPF];
  
  // Separate the quadrature states from the rest of the states
  MX yf_aug = xf_aug[range((ns+1)*ny_)];
  MX qf_aug = xf_aug[range((ns+1)*ny_,(ns+1)*nx_)];
  MX ypf_aug = xpf_aug[range((ns+1)*ny_)];
  MX qpf_aug = xpf_aug[range((ns+1)*ny_,(ns+1)*nx_)];
  
  // Get the state and state derivative at the final time
  MX yf = yf_aug[range(ny_)];
  MX qf = qf_aug[range(nq_)];
  MX xf = vertcat(yf,qf);
  MX ypf = yf_aug[range(ny_)];
  MX qpf = qf_aug[range(nq_)];
  MX xpf = vertcat(ypf,qpf);
  
  // Get the sensitivitiy equations state at the final time
  MX yf_sens = yf_aug[range(ny_,(ns+1)*ny_)];
  MX qf_sens = qf_aug[range(nq_,(ns+1)*nq_)];
  MX ypf_sens = yf_aug[range(ny_,(ns+1)*ny_)];
  MX qpf_sens = qf_aug[range(nq_,(ns+1)*nq_)];

  // Reshape the sensitivity state and state derivatives
  yf_sens = trans(reshape(yf_sens,ns,ny_));
  ypf_sens = trans(reshape(ypf_sens,ns,ny_));
  qf_sens = trans(reshape(qf_sens,ns,nq_));
  qpf_sens = trans(reshape(qpf_sens,ns,nq_));
  
  // Finally, we are able to get the Jacobian
  MX xf_sens = vertcat(yf_sens,qf_sens);
  MX xpf_sens = vertcat(ypf_sens,qpf_sens);

  MXFunction intjac(jac_in,xf_sens);
  
  return intjac;
  } else { 
    Integrator I;
    I.assignNode(this);
    return Jacobian(I,iind,oind);
  }
}

CRSSparsity IntegratorInternal::getJacSparsity(int iind, int oind){
  if(iind==INTEGRATOR_XP0){
    // Function value does not depend on the state derivative initial guess
    return CRSSparsity();
  } else {
    // Default (dense) sparsity
    return FXInternal::getJacSparsity(iind,oind);
  }
}

void IntegratorInternal::setInitialTime(double t0){
  t0_ = t0;
}

void IntegratorInternal::setFinalTime(double tf){
  tf_ = tf;
}


} // namespace CasADi


