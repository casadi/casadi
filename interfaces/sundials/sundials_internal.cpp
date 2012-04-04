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

#include "sundials_internal.hpp"

using namespace std;
namespace CasADi{
namespace Sundials{

SundialsInternal::SundialsInternal(const FX& f, const FX& q) : IntegratorInternal(f,q){
  addOption("max_num_steps",               OT_INTEGER, 10000); // maximum number of steps
  addOption("reltol",                      OT_REAL,    1e-6); // relative tolerence for the IVP solution
  addOption("abstol",                      OT_REAL,    1e-8); // absolute tolerence  for the IVP solution
  addOption("exact_jacobian",              OT_BOOLEAN,  false);
  addOption("upper_bandwidth",             OT_INTEGER); // upper band-width of banded jacobians
  addOption("lower_bandwidth",             OT_INTEGER); // lower band-width of banded jacobians
  addOption("linear_solver",               OT_STRING, "dense","","user_defined|dense|banded|iterative");
  addOption("iterative_solver",            OT_STRING, "gmres","","gmres|bcgstab|tfqmr");
  addOption("pretype",                     OT_STRING, "none","","none|left|right|both");
  addOption("max_krylov",                  OT_INTEGER,  10);        // maximum krylov subspace size
  addOption("is_differential",             OT_INTEGERVECTOR, GenericType(), "A vector with a boolean describing the nature for each state.");
  addOption("sensitivity_method",          OT_STRING,  "simultaneous","","simultaneous|staggered");
  addOption("max_multistep_order",         OT_INTEGER, 5);
  addOption("use_preconditioner",          OT_BOOLEAN, false); // precondition an iterative solver

  // Quadratures
  addOption("quad_err_con",                OT_BOOLEAN,false); // should the quadratures affect the step size control

  // Forward sensitivity problem
  addOption("fsens_err_con",               OT_BOOLEAN,true); // include the forward sensitivities in all error controls
  addOption("finite_difference_fsens",     OT_BOOLEAN, false); // use finite differences to approximate the forward sensitivity equations (if AD is not available)
  addOption("fsens_reltol",                OT_REAL); // relative tolerence for the forward sensitivity solution [default: equal to reltol]
  addOption("fsens_abstol",                OT_REAL); // absolute tolerence for the forward sensitivity solution [default: equal to abstol]
  addOption("fsens_scaling_factors",       OT_REALVECTOR); // scaling factor for the components if finite differences is used
  addOption("fsens_sensitiviy_parameters", OT_INTEGERVECTOR); // specifies which components will be used when estimating the sensitivity equations

  // Adjoint sensivity problem
  addOption("steps_per_checkpoint",        OT_INTEGER,20); // number of steps between two consecutive checkpoints
  addOption("interpolation_type",          OT_STRING,"hermite","type of interpolation for the adjoint sensitivities","hermite|polynomial");
  addOption("asens_upper_bandwidth",       OT_INTEGER); // upper band-width of banded jacobians
  addOption("asens_lower_bandwidth",       OT_INTEGER); // lower band-width of banded jacobians
  addOption("asens_linear_solver",         OT_STRING, "dense","","dense|banded|iterative");
  addOption("asens_iterative_solver",      OT_STRING, "gmres","","gmres|bcgstab|tfqmr");
  addOption("asens_pretype",               OT_STRING, "none","","none|left|right|both");
  addOption("asens_max_krylov",            OT_INTEGER,  10);        // maximum krylov subspace size
  addOption("asens_reltol",                OT_REAL); // relative tolerence for the adjoint sensitivity solution [default: equal to reltol]
  addOption("asens_abstol",                OT_REAL); // absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]
  addOption("linear_solver_creator",       OT_LINEARSOLVER, GenericType(), "An linear solver creator function");
  addOption("linear_solver_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the linear solver");
}

SundialsInternal::~SundialsInternal(){ 
}

void SundialsInternal::init(){
  // Call the base class method
  IntegratorInternal::init();
  
  // Read options
  abstol_ = getOption("abstol");
  reltol_ = getOption("reltol");
  exact_jacobian_ = getOption("exact_jacobian");
  max_num_steps_ = getOption("max_num_steps");
  finite_difference_fsens_ = getOption("finite_difference_fsens");
  fsens_abstol_ = hasSetOption("fsens_abstol") ? double(getOption("fsens_abstol")) : abstol_;
  fsens_reltol_ = hasSetOption("fsens_reltol") ? double(getOption("fsens_reltol")) : reltol_;
  asens_abstol_ = hasSetOption("asens_abstol") ? double(getOption("asens_abstol")) : abstol_;
  asens_reltol_ = hasSetOption("asens_reltol") ? double(getOption("asens_reltol")) : reltol_;
  stop_at_end_ = getOption("stop_at_end");
}

void SundialsInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
  IntegratorInternal::deepCopyMembers(already_copied);
}

} // namespace Sundials
} // namespace CasADi


