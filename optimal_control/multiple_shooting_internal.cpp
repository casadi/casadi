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
#include "../casadi/fx/sx_function.hpp"
#include "../casadi/matrix/matrix_tools.hpp"
#include "../casadi/sx/sx_tools.hpp"
#include "../casadi/mx/mx_tools.hpp"
#include "../casadi/stl_vector_tools.hpp"

namespace CasADi{
  namespace OptimalControl{
    
    
MultipleShootingInternal::MultipleShootingInternal(const FX& ffcn, const FX& mfcn, const FX& cfcn, const FX& rfcn) : OCPSolverInternal(ffcn, mfcn, cfcn, rfcn){
}

void MultipleShootingInternal::init(){
  // Initialize the base classes
  OCPSolverInternal::init();

  // Get final time
  double tf = getOption("final_time").toDouble();

  // Set time grid
  for(int k=0; k<=nk_; ++k)
    input(OCP_T).at(k) = (k*tf)/nk_;

  pF_ = Parallelizer(vector<FX>(nk_,ffcn_));
  pF_.setOption("save_corrected_input",true);
  pF_.init();
  
  Jacobian IjacX(ffcn_,INTEGRATOR_X0,INTEGRATOR_XF);
  pJX_ = Parallelizer(vector<FX>(nk_,IjacX));
  pJX_.setOption("save_corrected_input",true);
  pJX_.init();
  
  Jacobian IjacP(ffcn_,INTEGRATOR_P,INTEGRATOR_XF);
  pJP_ = Parallelizer(vector<FX>(nk_,IjacP));
  pJP_.setOption("save_corrected_input",true);
  pJP_.init();

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
  vector<MX> pI_in(nk_*INTEGRATOR_NUM_IN);
  for(int k=0; k<nk_; ++k){
    pI_in[INTEGRATOR_T0 + k*INTEGRATOR_NUM_IN] = input(OCP_T)[0]; // should be k
    pI_in[INTEGRATOR_TF + k*INTEGRATOR_NUM_IN] = input(OCP_T)[1]; // should be k+1
    pI_in[INTEGRATOR_P + k*INTEGRATOR_NUM_IN] = vertcat(P,U[k]);
    pI_in[INTEGRATOR_X0 + k*INTEGRATOR_NUM_IN] = X[k];
    pI_in[INTEGRATOR_Z0 + k*INTEGRATOR_NUM_IN] = Z[k];
    pI_in[INTEGRATOR_XP0+ k*INTEGRATOR_NUM_IN] = XP[k];
  }

  // Evaluate in parallel
  vector<MX> pI_out = pF_.call(pI_in);

  //Constraint function
  vector<MX> g(nk_);

  // Collect the outputs
  for(int k=0; k<nk_; ++k){
    // Get the output
    MX XF = pI_out[INTEGRATOR_XF + INTEGRATOR_NUM_OUT*k];
    
    //append continuity constraints
    g[k] = XF-X[k+1];
  }

  // df_dx blocks
  std::vector<Matrix<SX> > Jfx = symbolic("Jfx", nx_, nx_, nk_);
  
  // df_dp blocks
  std::vector<Matrix<SX> > Jfp = symbolic("Jfp", nx_, nu_, nk_);
  
  //Jacobian mapping
  SXMatrix Jgv(nx_*nk_,NV);
  for(int k=0; k<nk_; ++k){
    
    // Jacobian mapping, row by row
    for(int ii=0; ii<nx_; ++ii){
      // Row and column offset
      int i = nx_*k+ii;
      int j = (nx_+nu_)*k;
      
      // Add df/dxk
      for(int jj=0; jj<nx_; ++jj)
        Jgv(i, j+jj) = Jfx[k](ii,jj);
      
      // Add df_du
      for(int jj=0; jj<nu_; ++jj)
        Jgv(i, j+nx_+jj) = Jfp[k](ii,jj);
      
      // Add df/dx_{k+1}
      Jgv(i, j+nx_+nu_+ii) = -1;
    }
  }
  
  // All blocks
  vector<SXMatrix> Jall;
  Jall.insert(Jall.end(),Jfx.begin(),Jfx.end());
  Jall.insert(Jall.end(),Jfp.begin(),Jfp.end());
  
  J_mapping_ = SXFunction(Jall,Jgv);
  J_mapping_.init();

  // Create Jacobian function
  J_ = CFunction(MultipleShootingInternal::jacobian_wrapper);
  J_.setUserData(this);
  
  J_.setNumInputs(1);
  J_.input(0) = Matrix<double>(NV,1,0);

  J_.setNumOutputs(1);
  J_.output(0) = J_mapping_.output(0);
  
  // Objective function
  F_ = MXFunction(V,mfcn_.call(X[nk_-1]));

  // Terminal constraints
  G_ = MXFunction(V,vertcat(g));
}

void MultipleShootingInternal::constraint_wrapper(CFunction &f, int fsenk_order, int asenk_order, void* user_data){
  casadi_assert(fsenk_order==0 && asenk_order==0);
  MultipleShootingInternal* this_ = (MultipleShootingInternal*)user_data;
  this_->constraint(f,fsenk_order,asenk_order);
}

void MultipleShootingInternal::constraint(CFunction &f, int fsenk_order, int asenk_order){
  
}

void MultipleShootingInternal::jacobian_wrapper(CFunction &f, int fsenk_order, int asenk_order, void* user_data){
  casadi_assert(fsenk_order==0 && asenk_order==0);
  MultipleShootingInternal* this_ = (MultipleShootingInternal*)user_data;
  this_->jacobian(f,fsenk_order,asenk_order);
}

void MultipleShootingInternal::jacobian(CFunction &f, int fsenk_order, int asenk_order){
  // Functions that can be evaluated in parallel
  Parallelizer J[2] = {pJX_,pJP_};
  
  for(int j=0; j<2; ++j){
    // Pass inputs to parallelizer
    for(int k=0; k<nk_; ++k){
      J[j].setInput(input(OCP_T)[0], INTEGRATOR_T0 + k*INTEGRATOR_NUM_IN); // should be k
      J[j].setInput(input(OCP_T)[1], INTEGRATOR_TF + k*INTEGRATOR_NUM_IN); // should be k+1
      J[j].setInput(&f.input()[k*(nx_+nu_)], INTEGRATOR_X0 + k*INTEGRATOR_NUM_IN);
      J[j].setInput(&f.input()[k*(nx_+nu_)+nx_], INTEGRATOR_P  + k*INTEGRATOR_NUM_IN);
    }
    
    // Evaluate function
    J[j].evaluate();
    
    // Save to mapping
    for(int k=0; k<nk_; ++k)
      J_mapping_.setInput(J[j].output(k),j*nk_+k);
  }
  
  // Construct the full Jacobian
  J_mapping_.evaluate();
  for(int i=0; i<f.getNumOutputs(); ++i)
    J_mapping_.getOutput(f.output(i),i);
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

