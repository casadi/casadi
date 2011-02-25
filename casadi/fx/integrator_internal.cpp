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
#include "integrator_jacobian_internal.hpp"

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
  addOption("t0",                          OT_REAL, 0.0); // start of the integration (gives an initial value for INTEGRATOR_T0, which will be removed)
  addOption("tf",                          OT_REAL, 1.0); // end of the integration (gives an initial value for INTEGRATOR_TF, which will be removed)

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
  nz_ = 0;
}

IntegratorInternal::~IntegratorInternal(){ 
}

void IntegratorInternal::setDimensions(int nx, int np, int nz){
  nx_ = nx;
  np_ = np;
  nz_ = nz;
  
  // Allocate space for inputs
  input_.resize(INTEGRATOR_NUM_IN);
  input(INTEGRATOR_T0)  = DMatrix(1,1,0); // initial time
  input(INTEGRATOR_TF)  = DMatrix(1,1,0); // final time
  input(INTEGRATOR_X0)  = DMatrix(nx_,1,0); // initial state value
  input(INTEGRATOR_XP0) = DMatrix(nx_,1,0); // initial state derivative value
  input(INTEGRATOR_P)   = DMatrix(np_,1,0); // parameter
  input(INTEGRATOR_Z0)  = DMatrix(nz_,1,0); // initial algebraic statee
  
  // Allocate space for outputs
  output_.resize(INTEGRATOR_NUM_OUT);
  output(INTEGRATOR_XF) = DMatrix(nx_,1,0);
  output(INTEGRATOR_XPF)= DMatrix(nx_,1,0);
  output(INTEGRATOR_ZF) = DMatrix(nz_,1,0);
}

void IntegratorInternal::evaluate(int fsens_order, int asens_order){
  double t0 = input(INTEGRATOR_T0)[0];
  double tf = input(INTEGRATOR_TF)[0];
  
  // Reset solver
  reset(fsens_order, asens_order);

  // Set the stop time of the integration -- don't integrate past this point
  if(stop_at_end_)
    setStopTime(tf);

  // Advance solution in time
  integrate(tf);

  if(asens_order>0){
    // Re-initialize backward problem
    resetAdj();

    // Integrate backwards to the beginning
    integrateAdj(t0);
  }
}

void IntegratorInternal::init(){
  // Call the base class method
  FXInternal::init();

  // read options
  exact_jacobian_ = getOption("exact_jacobian").toInt();
  abstol_ = getOption("abstol").toDouble(); // TODO: change to vector tolerences
  reltol_ = getOption("reltol").toDouble();
  max_num_steps_ = getOption("max_num_steps").toInt();
  finite_difference_fsens_ = getOption("finite_difference_fsens").toInt();
  fsens_abstol_ = hasSetOption("fsens_abstol") ? getOption("fsens_abstol").toDouble() : abstol_;
  fsens_reltol_ = hasSetOption("fsens_reltol") ? getOption("fsens_reltol").toDouble() : reltol_;
  asens_abstol_ = hasSetOption("asens_abstol") ? getOption("asens_abstol").toDouble() : abstol_;
  asens_reltol_ = hasSetOption("asens_reltol") ? getOption("asens_reltol").toDouble() : reltol_;
  stop_at_end_ = getOption("stop_at_end").toInt();
  nrhs_ = getOption("nrhs").toInt();
  
  // Give an intial value for the time horizon
  input(INTEGRATOR_T0).set(getOption("t0").toDouble());
  input(INTEGRATOR_TF).set(getOption("tf").toDouble());
}

FX IntegratorInternal::jacobian(int iind, int oind){
  if(symbjac()){
  
  // Create a new integrator for the forward sensitivity equations
  Integrator fwdint = jac(iind,oind);

  // Number of sensitivities
  int ns;
  switch(iind){
    case INTEGRATOR_P: ns = np_; break;
    case INTEGRATOR_X0: ns = nx_; break;
    case INTEGRATOR_Z0: ns = nz_; break;
    default: casadi_assert_message(false,"iind must be INTEGRATOR_P, INTEGRATOR_X0 or INTEGRATOR_Z0");
  }
  
  // Generate an jacobian
  IntegratorJacobian intjac(fwdint);
  intjac->jacmap_ = jacmap(ns);
  
  // Derivative with respect to what?
  intjac.setOption("derivative_index",iind);
  
  return intjac;
  } else { 
    Integrator I;
    I.assignNode(this);
    return Jacobian(I,iind,oind);
  }
}



} // namespace CasADi


