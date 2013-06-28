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

#include "direct_single_shooting_internal.hpp"
#include "../symbolic/fx/integrator.hpp"
#include "../symbolic/matrix/matrix_tools.hpp"
#include "../symbolic/mx/mx_tools.hpp"
#include "../symbolic/stl_vector_tools.hpp"
#include "../symbolic/fx/fx_tools.hpp"

using namespace std;
namespace CasADi{
    
DirectSingleShootingInternal::DirectSingleShootingInternal(const FX& ffcn, const FX& mfcn, const FX& cfcn, const FX& rfcn) : OCPSolverInternal(ffcn, mfcn, cfcn, rfcn){
  addOption("parallelization", OT_STRING, GenericType(), "Passed on to CasADi::Parallelizer");
  addOption("nlp_solver",               OT_NLPSOLVER,  GenericType(), "An NLPSolver creator function");
  addOption("nlp_solver_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the NLP Solver");
  addOption("integrator",               OT_INTEGRATOR, GenericType(), "An integrator creator function");
  addOption("integrator_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the integrator");
}

DirectSingleShootingInternal::~DirectSingleShootingInternal(){
}

void DirectSingleShootingInternal::init(){
  // Initialize the base classes
  OCPSolverInternal::init();

  // Create an integrator instance
  integratorCreator integrator_creator = getOption("integrator");
  integrator_ = integrator_creator(ffcn_,FX());
  if(hasSetOption("integrator_options")){
    integrator_.setOption(getOption("integrator_options"));
  }

  // Set t0 and tf
  integrator_.setOption("t0",0);
  integrator_.setOption("tf",tf_/nk_);
  integrator_.init();
  
  // Path constraints present?
  bool path_constraints = nh_>0;
  
  // Count the total number of NLP variables
  int NV = np_ + // global parameters
           nx_ + // initial state
           nu_*nk_; // local control
           
  // Declare variable vector for the NLP
  // The structure is as follows:
  // np x 1  (parameters)
  // ------
  // nx x 1  (states at time i=0)
  // ------
  // nu x 1  (controls in interval i=0)
  // .....
  // nx x 1  (controls in interval i=nk-1)
  
  MX V = msym("V",NV);
  int offset = 0;

  // Global parameters
  MX P = V[Slice(0,np_)];
  offset += np_;

  // Initial state
  MX X0 = V[Slice(offset,offset+nx_)];
  offset += nx_;
  
  // Control for each shooting interval
  vector<MX> U(nk_);
  for(int k=0; k<nk_; ++k){ // interior nodes
    U[k] = V[range(offset,offset+nu_)];
    offset += nu_;
  }
  
  // Make sure that the size of the variable vector is consistent with the number of variables that we have referenced
  casadi_assert(offset==NV);

  // Current state
  MX X = X0;

  // Objective
  MX nlp_j = 0;

  // Constraints
  vector<MX> nlp_g;
  nlp_g.reserve(nk_*(path_constraints ? 2 : 1));
  
  // For all shooting nodes
  for(int k=0; k<nk_; ++k){
    // Integrate
    vector<MX> int_out = integrator_.call(integratorIn("x0",X,"p",vertcat(P,U[k])));

    // Store expression for state trajectory
    X = int_out[INTEGRATOR_XF];
    
    // Add constraints on the state
    nlp_g.push_back(X);

    // Add path constraints
    if(path_constraints){
      vector<MX> cfcn_out = cfcn_.call(daeIn("x",X,"p",U[k])); // TODO: Change signature of cfcn_: remove algebraic variable, add control
      nlp_g.push_back(cfcn_out.at(0));
    }
  }

  // Terminal cost
  MX jk = mfcn_.call(mayerIn("x",X,"p",P)).at(0);
  nlp_j += jk;

  // NLP
  nlp_ = MXFunction(nlpIn("x",V),nlpOut("f",nlp_j,"g",vertcat(nlp_g)));
  nlp_.setOption("name","nlp");
  nlp_.init();
    
  // Get the NLP creator function
  NLPSolverCreator nlp_solver_creator = getOption("nlp_solver");
  
  // Allocate an NLP solver
  nlp_solver_ = nlp_solver_creator(nlp_);
  
  // Pass options
  if(hasSetOption("nlp_solver_options")){
    const Dictionary& nlp_solver_options = getOption("nlp_solver_options");
    nlp_solver_.setOption(nlp_solver_options);
  }
  
  // Initialize the solver
  nlp_solver_.init();
}

void DirectSingleShootingInternal::getGuess(vector<double>& V_init) const{
  // OCP solution guess
  const Matrix<double> &p_init = input(OCP_P_INIT);
  const Matrix<double> &x_init = input(OCP_X_INIT);
  const Matrix<double> &u_init = input(OCP_U_INIT);
  
  // Running index
  int el=0;
  
  // Pass guess for parameters
  for(int i=0; i<np_; ++i){
    V_init[el++] = p_init.elem(i);
  }
  
  // Pass guess for the initial state
  for(int i=0; i<nx_; ++i){
    V_init[el++] = x_init.elem(i,0);
  }
  
  // Pass guess for control
  for(int k=0; k<nk_; ++k){
    for(int i=0; i<nu_; ++i){
      V_init[el++] = u_init.elem(i,k);
    }
  }
  
  casadi_assert(el==V_init.size());
}

void DirectSingleShootingInternal::getVariableBounds(vector<double>& V_min, vector<double>& V_max) const{
  // OCP variable bounds 
  const Matrix<double> &p_min = input(OCP_LBP);
  const Matrix<double> &p_max = input(OCP_UBP);
  const Matrix<double> &x_min = input(OCP_LBX);
  const Matrix<double> &x_max = input(OCP_UBX);
  const Matrix<double> &u_min = input(OCP_LBU);
  const Matrix<double> &u_max = input(OCP_UBU);


  // Running index
  int min_el=0, max_el=0;
  
  // Pass bounds on parameters
  for(int i=0; i<np_; ++i){
    V_min[min_el++] = p_min.elem(i);
    V_max[max_el++] = p_max.elem(i);
  }

  // Pass bounds on initial state
  for(int i=0; i<nx_; ++i){
    V_min[min_el++] = x_min.elem(i,0);
    V_max[max_el++] = x_max.elem(i,0);
  }
  
  // Pass bounds on control
  for(int k=0; k<nk_; ++k){
    for(int i=0; i<nu_; ++i){
      V_min[min_el++] = u_min.elem(i,k);
      V_max[max_el++] = u_max.elem(i,k);
    }
  }

  casadi_assert(min_el==V_min.size() && max_el==V_max.size());
}

void DirectSingleShootingInternal::getConstraintBounds(vector<double>& G_min, vector<double>& G_max) const{
  // OCP constraint bounds
  const Matrix<double> &x_min = input(OCP_LBX);
  const Matrix<double> &x_max = input(OCP_UBX);
  const Matrix<double> &h_min = input(OCP_LBH);
  const Matrix<double> &h_max = input(OCP_UBH);
  
  // Running index
  int min_el=0, max_el=0;
  
  for(int k=0; k<nk_; ++k){
    for(int i=0; i<nx_; ++i){
      G_min[min_el++] = x_min.elem(i,k+1);
      G_max[max_el++] = x_max.elem(i,k+1);
    }
    
    for(int i=0; i<nh_; ++i){
      G_min[min_el++] = h_min.elem(i,k);
      G_max[max_el++] = h_max.elem(i,k);
    }
  }
  casadi_assert(min_el==G_min.size() && max_el==G_max.size());
}

void DirectSingleShootingInternal::setOptimalSolution(const vector<double> &V_opt){
  // OCP solution
  Matrix<double> &p_opt = output(OCP_P_OPT);
  Matrix<double> &x_opt = output(OCP_X_OPT);
  Matrix<double> &u_opt = output(OCP_U_OPT);
  
  // Running index
  int el=0;

  // Pass optimized parameters
  for(int i=0; i<np_; ++i){
    p_opt(i) = V_opt[el++];
  }
    
  // Pass optimized initial state
  for(int i=0; i<nx_; ++i){
    x_opt(i,0) = V_opt[el++];
  }
  
  // Pass optimized control
  for(int k=0; k<nk_; ++k){
    for(int i=0; i<nu_; ++i){
      u_opt(i,k) = V_opt[el++];
    }
  }
  casadi_assert(el==V_opt.size());

  // Get the rest of the state trajectory
  const vector<double>& g_opt = nlp_solver_.output("g").data();

  // Loop over the constraints
  el = 0;
  for(int k=0; k<nk_; ++k){

    // Get the state trajectory
    for(int i=0; i<nx_; ++i){
      x_opt(i,k+1) = g_opt[el++];
    }
    
    // Skip the path constraints (for now)
    el += nh_;
  }
  casadi_assert(el==g_opt.size());
}

void DirectSingleShootingInternal::evaluate(int nfdir, int nadir){
  // get NLP variable bounds and initial guess
  getGuess(nlp_solver_.input(NLP_SOLVER_X0).data());
  getVariableBounds(nlp_solver_.input(NLP_SOLVER_LBX).data(),nlp_solver_.input(NLP_SOLVER_UBX).data());
       
  // get NLP constraint bounds
  getConstraintBounds(nlp_solver_.input(NLP_SOLVER_LBG).data(), nlp_solver_.input(NLP_SOLVER_UBG).data());
       
  //Solve the problem
  nlp_solver_.solve();
  
  // Save the optimal solution
  setOptimalSolution(nlp_solver_.output(NLP_SOLVER_X).data());

  // Save the optimal cost
  output(OCP_COST).set(nlp_solver_.output(NLP_SOLVER_F));
}


void DirectSingleShootingInternal::reportConstraints(std::ostream &stream) { 
  stream << "Reporting DirectSingleShooting constraints" << endl;
 
  CasADi::reportConstraints(stream,output(OCP_X_OPT),input(OCP_LBX),input(OCP_UBX), "states");
  CasADi::reportConstraints(stream,output(OCP_U_OPT),input(OCP_LBU),input(OCP_UBU), "controls");
  CasADi::reportConstraints(stream,output(OCP_P_OPT),input(OCP_LBP),input(OCP_UBP), "parameters");
 
}

} // namespace CasADi
