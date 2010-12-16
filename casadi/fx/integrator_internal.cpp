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

using namespace std;
namespace CasADi{

IntegratorInternal::IntegratorInternal(int nx, int np) : nx_(nx), np_(np){ 
  // set default options
  setOption("name","unnamed integrator"); // name of the function
  
  // Add new options
  addOption("max_num_steps",               OT_INTEGER, 10000); // maximum number of steps
  addOption("reltol",                      OT_REAL,    1e-6); // relative tolerence for the IVP solution
  addOption("abstol",                      OT_REAL,    1e-8); // absolute tolerence  for the IVP solution
  addOption("upper_bandwidth",             OT_INTEGER,  Option()); // upper band-width of banded jacobians
  addOption("lower_bandwidth",             OT_INTEGER,  Option()); // lower band-width of banded jacobians
  addOption("linear_solver",               OT_STRING, "dense"); // "dense", "banded" or "iterative"
  addOption("iterative_solver",            OT_STRING, "gmres"); // "gmres", "bcgstab", "tfqmr"
  addOption("pretype",                     OT_STRING, "none"); // "none", "left", "right", "both"
  addOption("exact_jacobian",              OT_BOOLEAN,  false);
  addOption("max_krylov",                  OT_INTEGER,  10);        // maximum krylov subspace size
  addOption("is_differential",             OT_INTEGERVECTOR,  Option());
  addOption("sensitivity_method",          OT_STRING,  "simultaneous"); // "simultaneous" or "staggered"
  addOption("max_multistep_order",         OT_INTEGER, 5);
  addOption("use_preconditioner",          OT_BOOLEAN, false); // precondition an iterative solver
  addOption("stop_at_end",                 OT_BOOLEAN, false); // Stop the integrator at the end of the interval

  // Quadratures
  addOption("quad_err_con",                OT_BOOLEAN,false); // should the quadratures affect the step size control

  // Forward sensitivity problem
  addOption("fsens_err_con",               OT_INTEGER, false); // include the forward sensitivities in all error controls
  addOption("finite_difference_fsens",     OT_BOOLEAN, false); // use finite differences to approximate the forward sensitivity equations (if AD is not available)
  addOption("fsens_reltol",                OT_REAL,    Option()); // relative tolerence for the forward sensitivity solution [default: equal to reltol]
  addOption("fsens_abstol",                OT_REAL,    Option()); // absolute tolerence for the forward sensitivity solution [default: equal to abstol]
  addOption("fsens_scaling_factors",       OT_REALVECTOR, Option()); // scaling factor for the components if finite differences is used
  addOption("fsens_sensitiviy_parameters", OT_INTEGERVECTOR, Option()); // specifies which components will be used when estimating the sensitivity equations

  // Adjoint sensivity problem
  addOption("steps_per_checkpoint",        OT_INTEGER,20); // number of steps between two consecutive checkpoints
  addOption("interpolation_type",          OT_STRING,"hermite"); // type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")
  addOption("asens_upper_bandwidth",       OT_INTEGER,  Option()); // upper band-width of banded jacobians
  addOption("asens_lower_bandwidth",       OT_INTEGER,  Option()); // lower band-width of banded jacobians
  addOption("asens_linear_solver",         OT_STRING, "dense"); // "dense", "banded" or "iterative"
  addOption("asens_iterative_solver",      OT_STRING, "gmres"); // "gmres", "bcgstab", "tfqmr"
  addOption("asens_pretype",               OT_STRING, "none"); // "none", "left", "right", "both"
  addOption("asens_max_krylov",            OT_INTEGER,  10);        // maximum krylov subspace size
  addOption("asens_reltol",                     OT_REAL,    Option()); // relative tolerence for the adjoint sensitivity solution [default: equal to reltol]
  addOption("asens_abstol",                     OT_REAL,    Option()); // absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]

  // Allocate space for inputs
  input_.resize(INTEGRATOR_NUM_IN);
  input_[INTEGRATOR_T0].setSize(1,1); // initial time
  input_[INTEGRATOR_TF].setSize(1,1); // final time
  input_[INTEGRATOR_X0].setSize(nx_,1); // initial state value
  input_[INTEGRATOR_XP0].setSize(nx_,1); // initial state derivative value
  input_[INTEGRATOR_P].setSize(np_,1); // parameter
  
  // Allocate space for outputs
  output_.resize(INTEGRATOR_NUM_OUT);
  output_[INTEGRATOR_XF].setSize(nx_,1);
  output_[INTEGRATOR_XPF].setSize(nx_,1);

}

IntegratorInternal::~IntegratorInternal(){ 
}

void IntegratorInternal::evaluate(int fsens_order, int asens_order){
  double t0 = input(INTEGRATOR_T0).data()[0];
  double tf = input(INTEGRATOR_TF).data()[0];
  
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
  FXNode::init();

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
}


} // namespace CasADi


