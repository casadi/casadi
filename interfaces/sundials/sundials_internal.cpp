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
#include "casadi/stl_vector_tools.hpp"
#include "casadi/matrix/matrix_tools.hpp"
#include "casadi/mx/mx_tools.hpp"
#include "casadi/sx/sx_tools.hpp"
#include "casadi/fx/mx_function.hpp"
#include "casadi/fx/sx_function.hpp"

INPUTSCHEME(IntegratorInput)
OUTPUTSCHEME(IntegratorOutput)

using namespace std;
namespace CasADi{
namespace Sundials{
  
SundialsInternal::SundialsInternal(const FX& f, const FX& g) : IntegratorInternal(f,g){
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
  addOption("sensitivity_method",          OT_STRING,  "simultaneous","","simultaneous|staggered");
  addOption("max_multistep_order",         OT_INTEGER, 5);
  addOption("use_preconditioner",          OT_BOOLEAN, false); // precondition an iterative solver
  addOption("stop_at_end",                 OT_BOOLEAN, false); // Stop the integrator at the end of the interval
  
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
  
  // Get the linear solver creator function
  if(linsol_.isNull() && hasSetOption("linear_solver_creator")){
    linearSolverCreator linear_solver_creator = getOption("linear_solver_creator");
  
    // Allocate an NLP solver
    linsol_ = linear_solver_creator(CRSSparsity());
  
    // Pass options
    if(hasSetOption("linear_solver_options")){
      const Dictionary& linear_solver_options = getOption("linear_solver_options");
      linsol_.setOption(linear_solver_options);
    }
  }
}

void SundialsInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
  IntegratorInternal::deepCopyMembers(already_copied);
  jac_ = deepcopy(jac_,already_copied);
  linsol_ = deepcopy(linsol_,already_copied);
}

SundialsIntegrator SundialsInternal::jac(bool with_x, bool with_p){
  // Make sure that we need sensitivities w.r.t. X0 or P (or both)
  casadi_assert(with_x || with_p);
  
  // Type cast to SXMatrix (currently supported)
  SXFunction f = shared_cast<SXFunction>(f_);
  if(f.isNull() != f_.isNull()) return SundialsIntegrator();
  
  // Number of state derivatives
  int n_xdot = f_.input(DAE_XDOT).numel();
  
  // Number of sensitivities
  int ns_x = with_x*nx_;
  int ns_p = with_p*np_;
  int ns = ns_x + ns_p;

  // Sensitivities and derivatives of sensitivities
  SXMatrix x_sens = ssym("x_sens",nx_,ns);
  SXMatrix z_sens = ssym("z_sens",nz_,ns);
  SXMatrix xdot_sens = ssym("xdot_sens",n_xdot,ns);
  
  // Directional derivative seeds
  vector<vector<SXMatrix> > fseed(ns);
  for(int d=0; d<ns; ++d){
    fseed[d].resize(DAE_NUM_IN);
    fseed[d][DAE_X] = x_sens(range(nx_),d);
    fseed[d][DAE_Z] = z_sens(range(nz_),d);
    fseed[d][DAE_P] = SXMatrix(f.inputSX(DAE_P).sparsity());
    if(with_p && d>=ns_x){
      fseed[d][DAE_P](d-ns_x) = 1;
    }
    fseed[d][DAE_T] = SXMatrix(f.inputSX(DAE_T).sparsity());
    if(n_xdot>0){
      fseed[d][DAE_XDOT] = xdot_sens(range(n_xdot),d);
    } else {
      fseed[d][DAE_XDOT] = SXMatrix(f.inputSX(DAE_XDOT).sparsity());
    }
  }
  
  // Calculate directional derivatives
  vector<vector<SXMatrix> > fsens(ns,f.outputsSX());
  vector<vector<SXMatrix> > aseedsens;
  f.evalSX(f.inputsSX(),const_cast<vector<SXMatrix>&>(f.outputsSX()),fseed,fsens,aseedsens,aseedsens,true);
  
  // Sensitivity equation (ODE)
  SXMatrix ode_aug = f.outputSX(DAE_ODE);
  for(int d=0; d<ns; ++d){
    ode_aug.append(fsens[d][DAE_ODE]);
  }

  // Sensitivity equation (ALG)
  SXMatrix alg_aug = f.outputSX(DAE_ALG);
  for(int d=0; d<ns; ++d){
    alg_aug.append(fsens[d][DAE_ALG]);
  }
  
  // Sensitivity equation (QUAD)
  SXMatrix quad_aug = f.outputSX(DAE_QUAD);
  for(int d=0; d<ns; ++d){
    quad_aug.append(fsens[d][DAE_QUAD]);
  }
  
  // Input arguments for the augmented DAE
  vector<SXMatrix> faug_in(DAE_NUM_IN);
  faug_in[DAE_T] = f.inputSX(DAE_T);
  faug_in[DAE_X] = vec(horzcat(f.inputSX(DAE_X),x_sens));
  if(nz_>0)    faug_in[DAE_Z] = vec(horzcat(f.inputSX(DAE_Z),z_sens));
  if(n_xdot>0) faug_in[DAE_XDOT] = vec(horzcat(f.inputSX(DAE_XDOT),xdot_sens));
  faug_in[DAE_P] = f.inputSX(DAE_P);
  
  // Create augmented DAE function
  SXFunction ffcn_aug(faug_in,daeOut(ode_aug,alg_aug,quad_aug));

  // Create integrator instance
  SundialsIntegrator integrator;
  integrator.assignNode(create(ffcn_aug,FX()));

  // Set options
  integrator.setOption(dictionary());
  integrator.setOption("nrhs",1+ns);
  
  // Pass linear solver
  if(!linsol_.isNull()){
    LinearSolver linsol_aug = shared_cast<LinearSolver>(linsol_.clone());
//    linsol_aug->nrhs_ = 1+ns; // FIXME!!!
//    integrator.setLinearSolver(linsol_aug,jac_); // FIXME!!!
    integrator.setLinearSolver(linsol_aug);
  }
  
  return integrator;
}

CRSSparsity SundialsInternal::getJacSparsity(int iind, int oind){
  // Default (dense) sparsity
    return FXInternal::getJacSparsity(iind,oind);
}

FX SundialsInternal::jacobian(const std::vector<std::pair<int,int> >& jblocks){
  bool with_x = false, with_p = false;
  for(int i=0; i<jblocks.size(); ++i){
    if(jblocks[i].second == INTEGRATOR_P){
      casadi_assert_message(jblocks[i].first == INTEGRATOR_XF,"IntegratorInternal::jacobian: Not derivative of state"); // can be removed?
      with_p = true;
    } else if(jblocks[i].second == INTEGRATOR_X0){
      casadi_assert_message(jblocks[i].first == INTEGRATOR_XF,"IntegratorInternal::jacobian: Not derivative of state"); // can be removed?
      with_x = true;
    }
  }
  
  // Create a new integrator for the forward sensitivity equations
  SundialsIntegrator fwdint = jac(with_x,with_p);
  
  // If creation was successful
  if(!fwdint.isNull()){
    fwdint.init();

    // Number of sensitivities
    int ns_x = with_x*nx_;
    int ns_p = with_p*np_;
    int ns = ns_x + ns_p;
    
    // Symbolic input of the Jacobian
    vector<MX> jac_in = symbolicInput();
    
    // Input to the augmented integrator
    vector<MX> fwdint_in(INTEGRATOR_NUM_IN);
    
    // Pass parameters without change
    fwdint_in[INTEGRATOR_P] = jac_in[INTEGRATOR_P];
    
    // Get the state
    MX x0 = jac_in[INTEGRATOR_X0];
    
    // Initial condition for the sensitivitiy equations
    DMatrix x0_sens(ns*nx_,1,0);
    
    if(with_x){
      for(int i=0; i<nx_; ++i){
	x0_sens.data()[i + i*ns_x] = 1;
      }
    }

    // Finally, we are ready to pass the initial condition for the state and state derivative
    fwdint_in[INTEGRATOR_X0] = vertcat(x0,MX(x0_sens));
    
    // Call the integrator with the constructed input (in fact, create a call node)
    vector<MX> fwdint_out = fwdint.call(fwdint_in);
    MX xf_aug = fwdint_out[INTEGRATOR_XF];
    MX qf_aug = fwdint_out[INTEGRATOR_QF];
    
    // Get the state and quadrature at the final time
    MX xf = xf_aug[range(nx_)];
    MX qf = qf_aug[range(nq_)];
    
    // Get the sensitivitiy equations state at the final time
    MX xf_sens = xf_aug[range(nx_,(ns+1)*nx_)];
    MX qf_sens = qf_aug[range(nq_,(ns+1)*nq_)];

    // Reshape the sensitivity state and state derivatives
    xf_sens = trans(reshape(xf_sens,ns,nx_));
    qf_sens = trans(reshape(qf_sens,ns,nq_));
    
    // Split up the Jacobians in parts for x0 and p
    MX J_xf_x0 = xf_sens(range(xf_sens.size1()),range(ns_x));
    MX J_xf_p  = xf_sens(range(xf_sens.size1()),range(ns_x,ns));
    MX J_qf_x0 = qf_sens(range(qf_sens.size1()),range(ns_x));
    MX J_qf_p  = qf_sens(range(qf_sens.size1()),range(ns_x,ns));
    
    // Output of the Jacobian
    vector<MX> jac_out(jblocks.size());
    for(int i=0; i<jblocks.size(); ++i){
      bool is_jac = jblocks[i].second >=0;
      bool is_x0 = jblocks[i].second==INTEGRATOR_X0;
      bool is_xf = jblocks[i].first==INTEGRATOR_XF;
      if(is_jac){
        if(is_x0){
          jac_out[i] = is_xf ? J_xf_x0 : J_qf_x0;
        } else {
          jac_out[i] = is_xf ? J_xf_p : J_qf_p;
        }
      } else {
        jac_out[i] = is_xf ? xf : qf;
      }
    }
    MXFunction intjac(jac_in,jac_out);

    return intjac;
  } else {
    return FXInternal::jacobian(jblocks);
  }
}

void SundialsInternal::setInitialTime(double t0){
  t0_ = t0;
}

void SundialsInternal::setFinalTime(double tf){
  tf_ = tf;
}

} // namespace Sundials
} // namespace CasADi


