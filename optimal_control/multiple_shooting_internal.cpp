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
  tf_ = getOption("final_time").toDouble();

  Jacobian IjacX(ffcn_,INTEGRATOR_X0,INTEGRATOR_XF);
  JX_ = Parallelizer(vector<FX>(nk_,IjacX));
  JX_.init();
  
  Jacobian IjacP(ffcn_,INTEGRATOR_P,INTEGRATOR_XF);
  JP_ = Parallelizer(vector<FX>(nk_,IjacP));
  JP_.init();

  //Number of discretized controls
  int NU = nu_*nk_;

  //Number of discretized states
  int NX = nx_*(nk_+1);

  //Declare variable vector
  int NV = NU+NX;
  MX V("V",NV);

  //Disretized control
  vector<MX> U;
  
  //Disretized state
  vector<MX> X;

  for(int k=0; k<nk_; ++k){
    X.push_back(V(range(k*(nx_+nu_),k*(nx_+nu_)+nx_),range(1)));
    U.push_back(V(range(k*(nx_+nu_)+nx_,k*(nx_+nu_)+nx_+nu_),range(1)));
  }
  X.push_back(V(range(nk_*(nx_+nu_),nk_*(nx_+nu_)+nx_),range(1)));
    
  //Beginning of each shooting interval (shift time horizon)
  vector<MX> T0(nk_,MX(0));

  //End of each shooting interval (shift time horizon)
  vector<MX> TF(nk_,MX(tf_/nk_));

  //The final state
  MX XF;
  
  //Constraint function with upper and lower bounds
  vector<MX> g;
  vector<SXMatrix> C;
  SXMatrix JJ(NX,NV);

  //Build up a graph of integrator calls
  for(int k=0; k<nk_; ++k){
    MX I_in[] = {T0[k],TF[k],X[k],U[k],MX(),MX()};
    
    //call the integrator
    XF = ffcn_.call(vector<MX>(I_in,I_in+6))[INTEGRATOR_XF];
    
    //append continuity constraints
    g.push_back(XF-X[k+1]);
    
    // Create a block
    stringstream ss;
    ss << k;
    
    SXMatrix CX = symbolic(string("X_")+ss.str(),nx_,nx_);
    SXMatrix CP = symbolic(string("P_")+ss.str(),nx_,nu_);
    C.push_back(CX);
    C.push_back(CP);
    for(int ii=0; ii<nx_; ++ii){
      for(int jj=0; jj<nx_; ++jj){
        JJ(nx_*k+ii, (nx_+nu_)*k+jj) = CX(ii,jj);
      }
      
      for(int jj=0; jj<nu_; ++jj){
        JJ(nx_*k+ii, (nx_+nu_)*k+nx_+jj) = CP(ii,jj);
      }
      
      JJ(nx_*k+ii, (nx_+nu_)*(k+1)+ii) = -1;
    }
  }
        
  J_mapping_ = SXFunction(C,JJ);
  J_mapping_.init();

  // Create Jacobian function
  J_ = CFunction(MultipleShootingInternal::jacobian_wrapper);
  J_.setUserData(this);
  
  J_.setNumInputs(1);
  J_.input(0) = Matrix<double>(NV,1,0);

  J_.setNumOutputs(1);
  J_.output(0) = J_mapping_.output(0);
  
  //State at the final time
  XF = X[nk_-1];

  //Objective function: L(T)
  F_ = MXFunction(V,mfcn_.call(XF));

  //Terminal constraints: 0<=[x(T);y(T)]<=0
  G_ = MXFunction(V,vertcat(g));
}
    
void MultipleShootingInternal::jacobian_wrapper(CFunction &f, int fsenk_order, int asenk_order, void* user_data){
  casadi_assert(fsenk_order==0 && asenk_order==0);
  MultipleShootingInternal* this_ = (MultipleShootingInternal*)user_data;
  this_->jacobian(f,fsenk_order,asenk_order);
}

void MultipleShootingInternal::jacobian(CFunction &f, int fsenk_order, int asenk_order){
  for(int i=0; i<nk_; ++i){
    JX_.input(INTEGRATOR_T0 + i*INTEGRATOR_NUM_IN)[0] = 0; 
    JX_.input(INTEGRATOR_TF + i*INTEGRATOR_NUM_IN)[0] = tf_/nk_; 
    for(int j=0; j<nx_; ++j)
      JX_.input(INTEGRATOR_X0 + i*INTEGRATOR_NUM_IN)[j] = f.input()[i*(nx_+nu_)+j];
    for(int j=0; j<nu_; ++j)
      JX_.input(INTEGRATOR_P  + i*INTEGRATOR_NUM_IN)[j] = f.input()[i*(nx_+nu_)+nx_+j];
  }

  for(int i=0; i<nk_; ++i){
    JP_.input(INTEGRATOR_T0 + i*INTEGRATOR_NUM_IN)[0] = 0; 
    JP_.input(INTEGRATOR_TF + i*INTEGRATOR_NUM_IN)[0] = tf_/nk_; 
    for(int j=0; j<nx_; ++j)
      JP_.input(INTEGRATOR_X0 + i*INTEGRATOR_NUM_IN)[j] = f.input()[i*(nx_+nu_)+j];
    for(int j=0; j<nu_; ++j)
      JP_.input(INTEGRATOR_P  + i*INTEGRATOR_NUM_IN)[j] = f.input()[i*(nx_+nu_)+nx_+j];
  }

  JX_.evaluate();
  JP_.evaluate();
  
  for(int i=0; i<nk_; ++i){
    J_mapping_.setInput(JX_.output(i),i*2);
    J_mapping_.setInput(JP_.output(i),i*2+1);
  }
  
  J_mapping_.evaluate();
  
  for(int i=0; i<f.getNumOutputs(); ++i){
    J_mapping_.getOutput(f.output(i),i);
  }
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
      V_init[k*(nu_+nx_)+i] = x_init[i+k*nx_];
      V_min[k*(nu_+nx_)+i] = x_min[i+k*nx_];
      V_max[k*(nu_+nx_)+i ]= x_max[i+k*nx_];
    }
  }
  
  // Pass bounds on control
  const vector<double> &u_min = input(OCP_LBU);
  const vector<double> &u_max = input(OCP_UBU);
  const vector<double> &u_init = input(OCP_U_INIT);
  for(int k=0; k<nk_; ++k){
    for(int i=0; i<nu_; ++i){
      V_min[k*(nu_+nx_)+nx_+i] = u_min[i+k*nu_];
      V_max[k*(nu_+nx_)+nx_+i] = u_max[i+k*nu_];
      V_init[k*(nu_+nx_)+nx_+i] = u_init[i+k*nu_];
    }
  }
  
  // Set constraint bounds to zero by default
  nlp_solver_.input(NLP_LBG).setZero();
  nlp_solver_.input(NLP_UBG).setZero();

  //Solve the problem
  nlp_solver_.solve();
}



  } // namespace OptimalControl
} // namespace CasADi

