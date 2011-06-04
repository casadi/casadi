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

#include "multiple_shooting_internal.hpp"
#include "../casadi/fx/integrator.hpp"
#include "../casadi/fx/jacobian.hpp"
#include "../casadi/matrix/matrix_tools.hpp"
#include "../casadi/mx/mx_tools.hpp"
#include "../casadi/stl_vector_tools.hpp"

using namespace std;
namespace CasADi{
  namespace OptimalControl{
    
    
MultipleShootingInternal::MultipleShootingInternal(const FX& ffcn, const FX& mfcn, const FX& cfcn, const FX& rfcn) : OCPSolverInternal(ffcn, mfcn, cfcn, rfcn){
  addOption("parallelization", OT_STRING);
}

MultipleShootingInternal::~MultipleShootingInternal(){
}

void MultipleShootingInternal::init(){
  // Initialize the base classes
  OCPSolverInternal::init();

  // Get final time
  double tf = getOption("final_time");

  // Path constraints present?
  bool path_constraints = nh_>0;
  
  // Set time grid
  for(int k=0; k<=nk_; ++k)
    input(OCP_T).at(k) = (k*tf)/nk_;

  // Count the total number of NLP variables
  int NV = np_ + // global parameters
           nx_*(nk_+1) + // local state
           nu_*nk_; // local control
           
  // Declare variable vector for the NLP
  MX V("V",NV);

  // Global parameters
  MX P = V(range(np_));

  // offset in the variable vector
  int v_offset=np_; 
  
  // Disretized variables for each shooting node
  vector<MX> X(nk_+1), U(nk_), XP(nk_);
  for(int k=0; k<=nk_; ++k){ // interior nodes
    // Local state
    X[k] = V[range(v_offset,v_offset+nx_)];
    v_offset += nx_;
    
    // Variables below do not appear at the end point
    if(k==nk_) break;
    
    // Local control
    U[k] = V[range(v_offset,v_offset+nu_)];
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
    int_in[k][INTEGRATOR_XP0] = XP[k];
  }

  // Input to the parallel function evaluation
  vector<vector<MX> > fcn_in(nk_);
  for(int k=0; k<nk_; ++k){
    fcn_in[k].resize(DAE_NUM_IN);
    fcn_in[k][DAE_T] = input(OCP_T)[k];
    fcn_in[k][DAE_P] = vertcat(P,U[k]);
    fcn_in[k][DAE_Y] = X[k];
    fcn_in[k][DAE_YDOT] = XP[k];
  }

  // Options for the parallelizer
  Dictionary paropt;
  paropt["save_corrected_input"] = true;
  
  // Transmit parallelization mode
  if(hasSetOption("parallelization"))
    paropt["parallelization"] = getOption("parallelization");
  
  // Evaluate function in parallel
  vector<vector<MX> > pI_out = ffcn_.call(int_in,paropt);

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
  G_ = MXFunction(V,g);

  // Objective function
  vector<MX> f = mfcn_.call(X.back());
  F_ = MXFunction(V,f);

  // Objective scaling factor
  MX sigma("sigma");
  
  // Lagrange multipliers
  MX lambda("lambda",g.size1());
  
  // Lagrangian
  MX L = sigma*f[0] + inner_prod(lambda,g);
  
  // Input of the function
  vector<MX> FG_in(3);
  FG_in[0] = V;
  FG_in[1] = lambda;
  FG_in[2] = sigma;

  // Output of the function
  vector<MX> FG_out(3);
  FG_out[0] = f[0];
  FG_out[1] = g;
  FG_out[2] = L;
  
  // Function that evaluates function and constraints that can also be used to get the gradient of the constraint
  FG_ = MXFunction(FG_in,FG_out);
  
  // Generate the Jacobian of the constraints
  G_.init();
  J_ = MXFunction(V,G_.jac());
}

void MultipleShootingInternal::getGuess(vector<double>& V_init) const{
  // OCP solution guess
  const Matrix<double> &p_init = input(OCP_P_INIT);
  const Matrix<double> &x_init = input(OCP_X_INIT);
  const Matrix<double> &u_init = input(OCP_U_INIT);
  
  // Running index
  int el=0;
  
    // Pass guess for parameters
    for(int i=0; i<np_; ++i){
      V_init[el++] = p_init(i);
    }
  
  for(int k=0; k<nk_; ++k){
    // Pass guess for state
    for(int i=0; i<nx_; ++i){
      V_init[el++] = x_init(i,k);
    }
    
    // Pass guess for control
    for(int i=0; i<nu_; ++i){
      V_init[el++] = u_init(i,k);
    }
  }

  // Pass guess for final state
  for(int i=0; i<nx_; ++i){
    V_init[el++] = x_init(i,nk_);
  }
  
  casadi_assert(el==V_init.size());
}

void MultipleShootingInternal::getVariableBounds(vector<double>& V_min, vector<double>& V_max) const{
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
    V_min[min_el++] = p_min(i);
    V_max[max_el++] = p_max(i);
  }

  for(int k=0; k<nk_; ++k){
    // Pass bounds on state
    for(int i=0; i<nx_; ++i){
      V_min[min_el++] = x_min(i,k);
      V_max[max_el++] = x_max(i,k);
    }
    
    // Pass bounds on control
    for(int i=0; i<nu_; ++i){
      V_min[min_el++] = u_min(i,k);
      V_max[max_el++] = u_max(i,k);
    }
  }

  // Pass bounds on final state
  for(int i=0; i<nx_; ++i){
    V_min[min_el++] = x_min(i,nk_);
    V_max[max_el++] = x_max(i,nk_);
  }
  
  std::cout << "42: " << min_el << std::endl;
  std::cout << "42: " << max_el << std::endl;
  std::cout << "42: " << V_min.size() << std::endl;
  std::cout << "42: " << V_max.size() << std::endl;
  casadi_assert(min_el==V_min.size() && max_el==V_max.size());
}

void MultipleShootingInternal::getConstraintBounds(vector<double>& G_min, vector<double>& G_max) const{
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
      G_min[min_el++] = h_min(i,k);
      G_max[max_el++] = h_max(i,k);
    }
  }
  casadi_assert(min_el==G_min.size() && max_el==G_max.size());
}

void MultipleShootingInternal::setOptimalSolution( const vector<double> &V_opt ){
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

void MultipleShootingInternal::evaluate(int nfdir, int nadir){
  // get NLP variable bounds and initial guess
  getGuess(nlp_solver_.input(NLP_X_INIT).data());
  getVariableBounds(nlp_solver_.input(NLP_LBX).data(),nlp_solver_.input(NLP_UBX).data());
       
  // get NLP constraint bounds
  getConstraintBounds(nlp_solver_.input(NLP_LBG).data(), nlp_solver_.input(NLP_UBG).data());
       
  //Solve the problem
  nlp_solver_.solve();
  
  // Save the optimal solution
  setOptimalSolution(nlp_solver_.output(NLP_X_OPT).data());
}



  } // namespace OptimalControl
} // namespace CasADi

