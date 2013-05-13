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

#include "direct_multiple_shooting_internal.hpp"
#include "../symbolic/fx/integrator.hpp"
#include "../symbolic/matrix/matrix_tools.hpp"
#include "../symbolic/mx/mx_tools.hpp"
#include "../symbolic/stl_vector_tools.hpp"
#include "../symbolic/fx/fx_tools.hpp"

using namespace std;
namespace CasADi{
    
DirectMultipleShootingInternal::DirectMultipleShootingInternal(const FX& ffcn, const FX& mfcn, const FX& cfcn, const FX& rfcn) : OCPSolverInternal(ffcn, mfcn, cfcn, rfcn){
  addOption("parallelization", OT_STRING, GenericType(), "Passed on to CasADi::Parallelizer");
  addOption("nlp_solver",               OT_NLPSOLVER,  GenericType(), "An NLPSolver creator function");
  addOption("nlp_solver_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the NLP Solver");
  addOption("integrator",               OT_INTEGRATOR, GenericType(), "An integrator creator function");
  addOption("integrator_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the integrator");
}

DirectMultipleShootingInternal::~DirectMultipleShootingInternal(){
}

void DirectMultipleShootingInternal::init(){
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
           nx_*(nk_+1) + // local state
           nu_*nk_; // local control
           
  // Declare variable vector for the NLP
  // The structure is as follows:
  // np x 1  (parameters)
  // ------
  // nx x 1  (states at time i=0)
  // nu x 1  (controls in interval i=0)
  // ------
  // nx x 1  (states at time i=1)
  // nu x 1  (controls in interval i=1)
  // ------
  // .....
  // ------
  // nx x 1  (states at time i=nk)
  
  MX V("V",NV);

  // Global parameters
  MX P = V(Slice(0,np_));

  // offset in the variable vector
  int v_offset=np_; 
  
  // Disretized variables for each shooting node
  vector<MX> X(nk_+1), U(nk_);
  for(int k=0; k<=nk_; ++k){ // interior nodes
    // Local state
    X[k] = V[Slice(v_offset,v_offset+nx_)];
    v_offset += nx_;
    
    // Variables below do not appear at the end point
    if(k==nk_) break;
    
    // Local control
    U[k] = V[Slice(v_offset,v_offset+nu_)];
    v_offset += nu_;
  }
  
  // Make sure that the size of the variable vector is consistent with the number of variables that we have referenced
  casadi_assert(v_offset==NV);

  // Input to the parallel integrator evaluation
  vector<vector<MX> > int_in(nk_);
  for(int k=0; k<nk_; ++k){
    int_in[k].resize(INTEGRATOR_NUM_IN);
    int_in[k][INTEGRATOR_P] = vertcat(P,U[k]);
    int_in[k][INTEGRATOR_X0] = X[k];
  }

  // Input to the parallel function evaluation
  vector<vector<MX> > fcn_in(nk_);
  for(int k=0; k<nk_; ++k){
    fcn_in[k].resize(DAE_NUM_IN);
    fcn_in[k][DAE_T] = (k*tf_)/nk_;
    fcn_in[k][DAE_P] = vertcat(P,U.at(k));
    fcn_in[k][DAE_X] = X[k];
  }

  // Options for the parallelizer
  Dictionary paropt;
  
  // Transmit parallelization mode
  if(hasSetOption("parallelization"))
    paropt["parallelization"] = getOption("parallelization");
  
  // Evaluate function in parallel
  vector<vector<MX> > pI_out = integrator_.call(int_in,paropt);

  // Evaluate path constraints in parallel
  vector<vector<MX> > pC_out;
  if(path_constraints)
    pC_out = cfcn_.call(fcn_in,paropt);
  
  //Constraint function
  vector<MX> gg(2*nk_);

  // Collect the outputs
  for(int k=0; k<nk_; ++k){
    //append continuity constraints
    gg[2*k] = pI_out[k][INTEGRATOR_XF] - X[k+1];
    
    // append the path constraints
    if(path_constraints)
      gg[2*k+1] = pC_out[k][0];
  }

  // Terminal constraints
  MX g = vertcat(gg);

  // Objective function
  MX f;
  if (mfcn_.getNumInputs()==1) {
    f = mfcn_.call(X.back()).front();
  } else {
    vector<MX> mfcn_argin(MAYER_NUM_IN); 
    mfcn_argin[MAYER_X] = X.back();
    mfcn_argin[MAYER_P] = P;
    f = mfcn_.call(mfcn_argin).front();
  }

  // NLP
  nlp_ = MXFunction(nlpIn("x",V),nlpOut("f",f,"g",g));
  nlp_.setOption("numeric_jacobian",false);
  nlp_.setOption("ad_mode","forward");
  nlp_.init();
  
  // Get the NLP creator function
  NLPSolverCreator nlp_solver_creator = getOption("nlp_solver");
  
  // Allocate an NLP solver
  nlp_solver_ = nlp_solver_creator(nlp_);
  
  // Symbolically calculate the gradient
  //nlp_solver_.setOption("generate_gradient",true);

  // Pass user options
  if(hasSetOption("nlp_solver_options")){
    const Dictionary& nlp_solver_options = getOption("nlp_solver_options");
    nlp_solver_.setOption(nlp_solver_options);
  }
  
  // Initialize the solver
  nlp_solver_.init();
}

void DirectMultipleShootingInternal::getGuess(vector<double>& V_init) const{
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
  
  for(int k=0; k<nk_; ++k){
    // Pass guess for state
    for(int i=0; i<nx_; ++i){
      V_init[el++] = x_init.elem(i,k);
    }
    
    // Pass guess for control
    for(int i=0; i<nu_; ++i){
      V_init[el++] = u_init.elem(i,k);
    }
  }
  
  // Pass guess for final state
  for(int i=0; i<nx_; ++i){
    V_init[el++] = x_init.elem(i,nk_);
  }
  
  casadi_assert(el==V_init.size());
}

void DirectMultipleShootingInternal::getVariableBounds(vector<double>& V_min, vector<double>& V_max) const{
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

  for(int k=0; k<nk_; ++k){
    // Pass bounds on state
    for(int i=0; i<nx_; ++i){
      V_min[min_el++] = x_min.elem(i,k);
      V_max[max_el++] = x_max.elem(i,k);
    }
    
    // Pass bounds on control
    for(int i=0; i<nu_; ++i){
      V_min[min_el++] = u_min.elem(i,k);
      V_max[max_el++] = u_max.elem(i,k);
    }
  }

  // Pass bounds on final state
  for(int i=0; i<nx_; ++i){
    V_min[min_el++] = x_min.elem(i,nk_);
    V_max[max_el++] = x_max.elem(i,nk_);
  }
  
  casadi_assert(min_el==V_min.size() && max_el==V_max.size());
}

void DirectMultipleShootingInternal::getConstraintBounds(vector<double>& G_min, vector<double>& G_max) const{
  // OCP constraint bounds
  const Matrix<double> &h_min = input(OCP_LBH);
  const Matrix<double> &h_max = input(OCP_UBH);
  
  // Running index
  int min_el=0, max_el=0;
  
  for(int k=0; k<nk_; ++k){
    for(int i=0; i<nx_; ++i){
      G_min[min_el++] = 0.;
      G_max[max_el++] = 0.;
    }
    
    for(int i=0; i<nh_; ++i){
      G_min[min_el++] = h_min.elem(i,k);
      G_max[max_el++] = h_max.elem(i,k);
    }
  }
  casadi_assert(min_el==G_min.size() && max_el==G_max.size());
}

void DirectMultipleShootingInternal::setOptimalSolution(const vector<double> &V_opt){
  // OCP solution
  Matrix<double> &p_opt = output(OCP_P_OPT);
  Matrix<double> &x_opt = output(OCP_X_OPT);
  Matrix<double> &u_opt = output(OCP_U_OPT);
  
  // Running index
  int el=0;

  // Pass optimized state
  for(int i=0; i<np_; ++i){
    p_opt(i) = V_opt[el++];
  }
    
  for(int k=0; k<nk_; ++k){
    
    // Pass optimized state
    for(int i=0; i<nx_; ++i){
      x_opt(i,k) = V_opt[el++];
    }
    
    // Pass optimized control
    for(int i=0; i<nu_; ++i){
      u_opt(i,k) = V_opt[el++];
    }
  }

  // Pass optimized terminal state
  for(int i=0; i<nx_; ++i){
    x_opt(i,nk_) = V_opt[el++];
  }
  casadi_assert(el==V_opt.size());
}

void DirectMultipleShootingInternal::evaluate(int nfdir, int nadir){
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


void DirectMultipleShootingInternal::reportConstraints(std::ostream &stream) { 
  stream << "Reporting DirectMultipleShooting constraints" << endl;
 
  CasADi::reportConstraints(stream,output(OCP_X_OPT),input(OCP_LBX),input(OCP_UBX), "states");
  CasADi::reportConstraints(stream,output(OCP_U_OPT),input(OCP_LBU),input(OCP_UBU), "controls");
  CasADi::reportConstraints(stream,output(OCP_P_OPT),input(OCP_LBP),input(OCP_UBP), "parameters");
 
}

} // namespace CasADi
