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

#include "integrator_jacobian_internal.hpp"
#include "integrator_internal.hpp"
#include "../stl_vector_tools.hpp"

using namespace std;
namespace CasADi{

IntegratorJacobianInternal::IntegratorJacobianInternal(const Integrator& integrator) : integrator_(integrator){
  addOption("derivative_index", OT_INTEGER, INTEGRATOR_P); // integrate with respect to what?
}

IntegratorJacobianInternal::~IntegratorJacobianInternal(){ 
}

void IntegratorJacobianInternal::init(){
  // Initialize the integrator
  integrator_.init();

  // Get the number of columns
  int nc = integrator_.getOption("nrhs").toInt();
  
  // Number of sensitivities
  ns_ = nc-1;
  
  // Number of states
  nx_ = integrator_->nx_/nc;

  // Number of parameters
  int np = integrator_->np_;

  // Set the dimensions
  input_.resize(INTEGRATOR_NUM_IN);
  input(INTEGRATOR_T0)     = DMatrix(1,1,0); // initial time
  input(INTEGRATOR_TF)     = DMatrix(1,1,0); // final time
  input(INTEGRATOR_X0)     = DMatrix(nx_,1,0); // initial state value
  input(INTEGRATOR_XP0)    = DMatrix(nx_,1,0); // initial state derivative value
  input(INTEGRATOR_P)      = DMatrix(np,1,0); // parameter
  
  // Allocate space for outputs
  output_.resize(1+INTEGRATOR_NUM_OUT);
  output(0)                = DMatrix(nx_,ns_,0);
  output(1+INTEGRATOR_XF)  = DMatrix(nx_,1,0);
  output(1+INTEGRATOR_XPF) = DMatrix(nx_,1,0);

  // Map Jacobian indices
  if(integrator_.hasSetOption("jacmap")){
    jacmap_ = integrator_.getOption("jacmap").toIntVector();
  } else {
    jacmap_.clear();
    jacmap_.resize(nx_*nc);
    for(int i=0; i<nx_*nc; ++i)
      jacmap_[i] = i;
  }

  // Initial value to the Jacobian
  if(integrator_.hasSetOption("jacinit")){
    jacinit_ = integrator_.getOption("jacinit").toDoubleVector();
  } else {
    jacinit_.clear();
    jacinit_.resize(nx_*ns_,0.0);
  }

  // Call the base class method
  FXInternal::init();
}

void IntegratorJacobianInternal::evaluate(int fsens_order, int asens_order){
  
  // Pass arguments to the integrator
  integrator_.setInput(input(INTEGRATOR_T0),INTEGRATOR_T0);
  integrator_.setInput(input(INTEGRATOR_TF),INTEGRATOR_TF);
  integrator_.setInput(input(INTEGRATOR_P),INTEGRATOR_P);

  // Initial value for the state
  const vector<double>& x0 = input(INTEGRATOR_X0);
  vector<double>& x0s = integrator_.input(INTEGRATOR_X0);
  fill(x0s.begin(),x0s.end(),0.0);
  for(int i=0; i<nx_; ++i)
    x0s[jacmap_[i]] = x0[i];

  // Initial values for the state derivatives
  for(int i=0; i<nx_; ++i)
    for(int j=0; j<ns_; ++j)
      x0s[jacmap_[nx_+j+i*ns_]] = jacinit_[j+i*ns_];
  
  // State derivative
  const vector<double>& xp0 = input(INTEGRATOR_XP0);
  vector<double>& xp0s = integrator_.input(INTEGRATOR_XP0);
  fill(xp0s.begin(),xp0s.end(),0.0);
  for(int i=0; i<nx_; ++i)
    xp0s[jacmap_[i]] = xp0[i];
  
  // Pass adjoint seeds
  if(asens_order>0){
    for(int dir=0; dir<nadir_; ++dir){
      const vector<double>& jacseed = adjSeed(0,dir);
      vector<double>& jacseed_s = integrator_.adjSeed(INTEGRATOR_XF,dir);
      for(int i=0; i<nx_; ++i)
        for(int j=0; j<ns_; ++j)
          jacseed_s[jacmap_[nx_+j+i*ns_]] = jacseed[j+i*ns_];
    }
  }
  
  // Solve the initial value problem
  integrator_.evaluate(fsens_order,asens_order);
  
  // Get the results
  const vector<double>& xfs = integrator_.output(INTEGRATOR_XF);
  const vector<double>& xpfs = integrator_.output(INTEGRATOR_XPF);

  // Jacobian
  vector<double>& jac = output(0);
  for(int i=0; i<nx_; ++i)
    for(int j=0; j<ns_; ++j)
      jac[j+i*ns_] = xfs[jacmap_[nx_+j+i*ns_]];
  
  // State
  vector<double>& xf = output(1+INTEGRATOR_XF);
  for(int i=0; i<nx_; ++i)
    xf[i] = xfs[jacmap_[i]];
    
  // State derivative
  vector<double>& xpf = output(1+INTEGRATOR_XPF);
  for(int i=0; i<nx_; ++i)
    xpf[i] = xpfs[jacmap_[i]];
  
  
  // Get adjoint sensitivities
  if(asens_order>0){
    for(int dir=0; dir<nadir_; ++dir){
      vector<double>& asens = adjSens(INTEGRATOR_X0,dir);
      fill(asens.begin(),asens.end(),0);
      
      const vector<double>& jacsens_s = integrator_.adjSens(INTEGRATOR_X0,dir);
      for(int i=0; i<nx_; ++i)
        for(int j=0; j<ns_; ++j)
          asens[i] += jacsens_s[jacmap_[nx_+j+i*ns_]];
        
      adjSens(INTEGRATOR_P,dir).set(integrator_.adjSens(INTEGRATOR_P,dir));
    }
  }
}

IntegratorJacobianInternal* IntegratorJacobianInternal::clone() const{
  // Return a deep copy
  IntegratorJacobianInternal* node = new IntegratorJacobianInternal(*this);
  node->integrator_ = shared_cast<Integrator>(integrator_.clone());
  return node;
}


} // namespace CasADi


