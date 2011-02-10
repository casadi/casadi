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

#include "fx_tools.hpp"
#include "integrator.hpp"
#include "jacobian.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../mx/mx_tools.hpp"

namespace CasADi{
    
MultipleShooting::MultipleShooting(const FX& fcn, int ns, int nx, int nu) : fcn_(fcn), ns_(ns), nx_(nx), nu_(nu){
  u_min_.resize(nu);   u_max_.resize(nu);   u_init_.resize(nu);
  x_min_.resize(nx);   x_max_.resize(nx);   x_init_.resize(nx);
  x0_min_.resize(nx);  x0_max_.resize(nx);
  xf_min_.resize(nx);  xf_max_.resize(nx);
}

void MultipleShooting::init(){
  Jacobian IjacX(fcn_,INTEGRATOR_X0,INTEGRATOR_XF);
  JX_ = Parallelizer(vector<FX>(ns_,IjacX));
  
  Jacobian IjacP(fcn_,INTEGRATOR_P,INTEGRATOR_XF);
  JP_ = Parallelizer(vector<FX>(ns_,IjacP));
  
  JX_.init();
  JP_.init();

  //Number of discretized controls
  int NU = nu_*ns_;

  //Number of discretized states
  int NX = nx_*(ns_+1);

  //Declare variable vector
  int NV = NU+NX;
  MX V("V",NV);

  //Disretized control
  vector<MX> U;
  
  //Disretized state
  vector<MX> X;

  for(int k=0; k<ns_; ++k){
    X.push_back(V(range(k*(nx_+nu_),k*(nx_+nu_)+nx_),range(1)));
    for(int i=0; i<nx_; ++i){
      V_init_.push_back(x_init_[i]);
    }
    
    if(k==0){
      for(int i=0; i<nx_; ++i){
        V_min_.push_back(x0_min_[i]);
        V_max_.push_back(x0_max_[i]);
      }
    } else {
      for(int i=0; i<nx_; ++i){
        V_min_.push_back(x_min_[i]);
        V_max_.push_back(x_max_[i]);
      }
    }
    
    U.push_back(V(range(k*(nx_+nu_)+nx_,k*(nx_+nu_)+nx_+nu_),range(1)));
    for(int i=0; i<nu_; ++i){
      V_min_.push_back(u_min_[i]);
      V_max_.push_back(u_max_[i]);
      V_init_.push_back(u_init_[i]);
    }
  }
  
  X.push_back(V(range(ns_*(nx_+nu_),ns_*(nx_+nu_)+nx_),range(1)));
  
  for(int i=0; i<nx_; ++i){
    V_init_.push_back(x_init_[i]);
    V_min_.push_back(x_min_[i]);
    V_max_.push_back(x_max_[i]);
  }
    
  //Beginning of each shooting interval (shift time horizon)
  vector<MX> T0(ns_,MX(0));

  //End of each shooting interval (shift time horizon)
  vector<MX> TF(ns_,MX(tf_/ns_));

  //The final state
  MX XF;
  
  //Constraint function with upper and lower bounds
  vector<MX> g;
  vector<SXMatrix> C;
  SXMatrix JJ(NX,NV);

  //Build up a graph of integrator calls
  for(int k=0; k<ns_; ++k){
    MX I_in[] = {T0[k],TF[k],X[k],U[k],MX(),MX()};
    
    //call the integrator
    XF = fcn_.call(vector<MX>(I_in,I_in+6))[INTEGRATOR_XF];
    
    //append continuity constraints
    g.push_back(XF-X[k+1]);
    for(int j=0; j<XF.numel(); ++j){
      G_min_.push_back(0);
      G_max_.push_back(0);
    }
    
    
    // Create a block
    stringstream ss;
    ss << k;
    
    SXMatrix CX = symbolic(string("X_")+ss.str(),nx_,nx_);
    SXMatrix CP = symbolic(string("P_")+ss.str(),nx_,nu_);
    C.push_back(CX);
    C.push_back(CP);
    for(int ii=0; ii<nx_; ++ii){
      for(int jj=0; jj<nx_; ++jj){
        JJ(nx_*k+ii, (nx_+nu_)*k+jj) = SX(CX(ii,jj)); // syntax workaround
      }
      
      for(int jj=0; jj<nu_; ++jj){
        JJ(nx_*k+ii, (nx_+nu_)*k+nx_+jj) = SX(CP(ii,jj)); // error!
      }
      
      JJ(nx_*k+ii, (nx_+nu_)*(k+1)+ii) = -1;
    }
  }
        
  J_mapping_ = SXFunction(C,JJ);
  J_mapping_.init();

  // Create Jacobian function
  J_ = CFunction(MultipleShooting::jacobian_wrapper);
  J_.setUserData(this);
  J_.setNumInputs(J_mapping_.getNumInputs());
  for(int i=0; i<J_.getNumInputs(); ++i){
    J_.input(i) = J_mapping_.input(i);
  }
  J_.setNumOutputs(J_mapping_.getNumOutputs());
  for(int i=0; i<J_.getNumOutputs(); ++i){
    J_.output(i) = J_mapping_.output(i);
  }
  
  //State at the final time
  XF = X[ns_-1];

  //Objective function: L(T)
  F_ = MXFunction(V,XF[2]);

  //Terminal constraints: 0<=[x(T);y(T)]<=0
  G_ = MXFunction(V,vertcat(g));
  
  // BUG: cannot remove
  Jacobian JJJ(G_); 
}
    
void MultipleShooting::jacobian_wrapper(CFunction &f, int fsens_order, int asens_order, void* user_data){
  casadi_assert(fsens_order==0 && asens_order==0);
  MultipleShooting* this_ = (MultipleShooting*)user_data;
  this_->jacobian(f,fsens_order,asens_order);
}

void MultipleShooting::jacobian(CFunction &f, int fsens_order, int asens_order){
  for(int i=0; i<ns_; ++i){
    JX_.input(INTEGRATOR_T0 + i*INTEGRATOR_NUM_IN)[0] = 0; 
    JX_.input(INTEGRATOR_TF + i*INTEGRATOR_NUM_IN)[0] = tf_/ns_; 
    for(int j=0; j<nx_; ++j)
      JX_.input(INTEGRATOR_X0 + i*INTEGRATOR_NUM_IN)[j] = f.input()[i*(nx_+nu_)+j];
    for(int j=0; j<nu_; ++j)
      JX_.input(INTEGRATOR_P  + i*INTEGRATOR_NUM_IN)[j] = f.input()[i*(nx_+nu_)+nx_+j];
  }

  for(int i=0; i<ns_; ++i){
    JP_.input(INTEGRATOR_T0 + i*INTEGRATOR_NUM_IN)[0] = 0; 
    JP_.input(INTEGRATOR_TF + i*INTEGRATOR_NUM_IN)[0] = tf_/ns_; 
    for(int j=0; j<nx_; ++j)
      JP_.input(INTEGRATOR_X0 + i*INTEGRATOR_NUM_IN)[j] = f.input()[i*(nx_+nu_)+j];
    for(int j=0; j<nu_; ++j)
      JP_.input(INTEGRATOR_P  + i*INTEGRATOR_NUM_IN)[j] = f.input()[i*(nx_+nu_)+nx_+j];
  }

  JX_.evaluate();
  JP_.evaluate();
  
  for(int i=0; i<ns_; ++i){
    J_mapping_.setInput(JX_.output(i),i*2);
    J_mapping_.setInput(JP_.output(i),i*2+1);
  }
  
  J_mapping_.evaluate();
  
  for(int i=0; i<f.getNumOutputs(); ++i){
    J_mapping_.getOutput(f.output(i),i);
  }
}
  
} // namespace CasADi

