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

namespace CasADi{
  namespace OptimalControl{
    
    
MultipleShootingInternal::MultipleShootingInternal(const FX& ffcn, const FX& mfcn, const FX& cfcn, const FX& rfcn) : OCPSolverInternal(ffcn, mfcn, cfcn, rfcn){
  addOption("parallelization", OT_STRING);
}

void MultipleShootingInternal::init(){
  // Initialize the base classes
  OCPSolverInternal::init();

  // Get final time
  double tf = getOption("final_time").toDouble();

  // Set time grid
  for(int k=0; k<=nk_; ++k)
    input(OCP_T).at(k) = (k*tf)/nk_;

  //Declare variable vector
  int NV = np_+nu_*nk_+nx_*(nk_+1);
  MX V("V",NV);
  
  //Disretized control
  vector<MX> U(nk_);
  for(int k=0; k<nk_; ++k)
    U[k] = V(range(np_+k*(nx_+nu_)+nx_,np_+k*(nx_+nu_)+nx_+nu_),range(1));
  
  //Disretized state
  vector<MX> X(nk_+1);
  for(int k=0; k<=nk_; ++k)
    X[k] = V(range(np_+k*(nx_+nu_),np_+k*(nx_+nu_)+nx_),range(1));

  // Algebraic state, state derivative
  vector<MX> Z(nk_), XP(nk_);

  // Parameters
  MX P = V(range(np_),range(1));

  // Input to the parallel evaluation
  vector<vector<MX> > ffcn_in(nk_);
  for(int k=0; k<nk_; ++k){
    ffcn_in[k].resize(INTEGRATOR_NUM_IN);
    ffcn_in[k][INTEGRATOR_T0] = input(OCP_T)[0]; // should be k
    ffcn_in[k][INTEGRATOR_TF] = input(OCP_T)[1]; // should be k+1
    ffcn_in[k][INTEGRATOR_P] = vertcat(P,U[k]);
    ffcn_in[k][INTEGRATOR_X0] = X[k];
    ffcn_in[k][INTEGRATOR_Z0] = Z[k];
    ffcn_in[k][INTEGRATOR_XP0] = XP[k];
  }

  // Options for the parallelizer
  Dictionary paropt;
  paropt["save_corrected_input"] = true;
  
  // Transmit parallelization mode
  if(hasSetOption("parallelization"))
    paropt["parallelization"] = getOption("parallelization");
  
  // Evaluate function in parallel
  vector<vector<MX> > pI_out = ffcn_.call(ffcn_in,paropt);

  //Constraint function
  vector<MX> g(nk_);

  // Collect the outputs
  for(int k=0; k<nk_; ++k){
    // Get the output
    MX XF = pI_out[k][INTEGRATOR_XF];
    
    //append continuity constraints
    g[k] = XF-X[k+1];
  }

  // Terminal constraints
  G_ = MXFunction(V,vertcat(g));

  // Objective function
  F_ = MXFunction(V,mfcn_.call(X[nk_-1]));

  // Evaluate Jacobian blocks in parallel
  Jacobian IjacX(ffcn_,INTEGRATOR_X0,INTEGRATOR_XF);
  vector<vector<MX> > pJX_out = IjacX.call(ffcn_in,paropt);
  
  Jacobian IjacP(ffcn_,INTEGRATOR_P,INTEGRATOR_XF);
  vector<vector<MX> > pJP_out = IjacP.call(ffcn_in,paropt);

  // Empty matrices with the same size as the Jacobian blocks
  MX e_JX = MX(pJX_out[0][0].size1(),pJX_out[0][0].size2());
  MX e_JP = MX(pJP_out[0][0].size1(),pJP_out[0][0].size2());
  
  // Identity matrices with the same size as the Jacobian blocks
  MX mi_JX = -eye<double>(IjacX.output().size1());
  
  // Vector of row blocks
  vector<MX> JJ;
  JJ.reserve(nk_);
  
  // Loop over the row blocks
  for(int i=0; i<nk_; ++i){
    
    // Block row of the Jacobian
    vector<MX> JJ_row;
    JJ_row.reserve(2*nk_+1);
    
    // Add all the blocks
    for(int j=0; j<nk_; ++j){
      if(j==i){                       // Diagonal block
        JJ_row.push_back(pJX_out[i][0]); // df_k/dx_k block
        JJ_row.push_back(pJP_out[i][0]); // df_k/du_k block
      } else if(j==i+1) {             // First upper subdiagonal block
        JJ_row.push_back(mi_JX);      // df_k/dx_{k+1} block
        JJ_row.push_back(e_JP);
      } else {                        // All other blocks
        JJ_row.push_back(e_JX);
        JJ_row.push_back(e_JP);
      }
    }
    
    if(i==nk_-1) {
      JJ_row.push_back(mi_JX);        // df_k/dx_{k+1} block
    } else {
      JJ_row.push_back(e_JX);
    }
    
    // Add to vector of rows
    JJ.push_back(horzcat(JJ_row));
  }
  
  // Create function
  J_ = MXFunction(V,vertcat(JJ));
}

void MultipleShootingInternal::evaluate(int fsenk_order, int asenk_order){
  // NLP Variable bounds and initial guess
  vector<double> &V_min = nlp_solver_.input(NLP_LBX);
  vector<double> &V_max = nlp_solver_.input(NLP_UBX);
  vector<double> &V_init = nlp_solver_.input(NLP_X_INIT);
  
  // Pass bounds on state
  const vector<double> &x_min = input(OCP_LBX);
  const vector<double> &x_max = input(OCP_UBX);
  const vector<double> &x_init = input(OCP_X_INIT);
  for(int k=0; k<nk_+1; ++k){
    for(int i=0; i<nx_; ++i){
      V_init[k*(nu_+nx_)+i] = x_init[k+i*(nk_+1)];
      V_min[k*(nu_+nx_)+i] = x_min[k+i*(nk_+1)];
      V_max[k*(nu_+nx_)+i ]= x_max[k+i*(nk_+1)];
    }
  }
  
  // Pass bounds on control
  const vector<double> &u_min = input(OCP_LBU);
  const vector<double> &u_max = input(OCP_UBU);
  const vector<double> &u_init = input(OCP_U_INIT);
  for(int k=0; k<nk_; ++k){
    for(int i=0; i<nu_; ++i){
      V_min[k*(nu_+nx_)+nx_+i] = u_min[k+i*nk_];
      V_max[k*(nu_+nx_)+nx_+i] = u_max[k+i*nk_];
      V_init[k*(nu_+nx_)+nx_+i] = u_init[k+i*nk_];
    }
  }
  
  // Set constraint bounds to zero by default
  nlp_solver_.input(NLP_LBG).setZero();
  nlp_solver_.input(NLP_UBG).setZero();

  //Solve the problem
  nlp_solver_.solve();

  // Optimal solution
  const vector<double> &V_opt = nlp_solver_.output(NLP_X_OPT);
  
  // Pass bounds on state
  vector<double> &x_opt = output(OCP_X_OPT);
  for(int k=0; k<nk_+1; ++k){
    for(int i=0; i<nx_; ++i){
      x_opt[k+i*(nk_+1)] = V_opt[k*(nu_+nx_)+i];
    }
  }
  
  // Pass bounds on control
  vector<double> &u_opt = output(OCP_U_OPT);
  for(int k=0; k<nk_; ++k){
    for(int i=0; i<nu_; ++i){
      u_opt[k+i*nk_] = V_opt[k*(nu_+nx_)+nx_+i];
    }
  }
}



  } // namespace OptimalControl
} // namespace CasADi

