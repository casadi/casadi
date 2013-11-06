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

#include "derivative_internal.hpp"
#include "../stl_vector_tools.hpp"
#include "../matrix/sparsity_tools.hpp"
#include "fx_internal.hpp"

using namespace std;

namespace CasADi{
  
DerivativeInternal::DerivativeInternal(const FX& fcn, int nfwd, int nadj) : fcn_(fcn), nfwd_(nfwd), nadj_(nadj){
    
  casadi_assert_message(nfwd==0 || int(fcn.getOption("max_number_of_fwd_dir")) > 0, "DerivativeInternal::DerivativeInternal: unable to request any forward sensitivities.");
  casadi_assert_message(nadj==0 || int(fcn.getOption("max_number_of_adj_dir")) > 0, "FXInternal::requestNumSens: unable to request any adjoint sensitivities.");
    
  casadi_assert(nfwd>0 || nadj>0);
}
  
DerivativeInternal::~DerivativeInternal(){
}

DerivativeInternal* DerivativeInternal::clone() const{
  return new DerivativeInternal(*this);
}

void DerivativeInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
  FXInternal::deepCopyMembers(already_copied);
  fcn_ = deepcopy(fcn_,already_copied);
}

void DerivativeInternal::init(){
  // Initialize the function if not already initialized
  if(!fcn_.isInit()) fcn_.init();

  // Number inputs and outputs of the function
  int num_in = fcn_.getNumInputs();
  int num_out = fcn_.getNumOutputs();
  
  // Allocate inputs
  setNumInputs(num_in*(1+nfwd_)+num_out*nadj_);
  int ind=0;
  for(int dir=-1; dir<nfwd_; ++dir){
    for(int i=0; i<num_in; ++i){
      input(ind++) = fcn_.input(i);
    }
  }
  for(int dir=0; dir<nadj_; ++dir){
    for(int i=0; i<num_out; ++i){
      input(ind++) = fcn_.output(i);
    }
  }
  casadi_assert(ind==getNumInputs());

  // Allocate outputs
  setNumOutputs(num_out*(1+nfwd_)+num_in*nadj_);
  ind = 0;
  for(int dir=-1; dir<nfwd_; ++dir){
    for(int i=0; i<num_out; ++i){
      output(ind++) = fcn_.output(i);
    }
  }
  for(int dir=0; dir<nadj_; ++dir){
    for(int i=0; i<num_in; ++i){
      output(ind++) = fcn_.input(i);
    }
  }
  casadi_assert(ind==getNumOutputs());
  
  // Call the base class init routine
  FXInternal::init();

  // Get the number of scalar inputs and outputs
  int num_in_scalar = fcn_.getNumScalarInputs();
  int num_out_scalar = fcn_.getNumScalarOutputs();
  
  // Adjoint mode penalty factor (adjoint mode is usually more expensive to calculate)
  int adj_penalty = 2;

  // Crude estimate of the cost of calculating the full Jacobian
  int full_jac_cost = std::min(num_in_scalar, adj_penalty*num_out_scalar);

  // Crude estimate of the cost of calculating the directional derivatives
  int der_dir_cost = nfwd_ + adj_penalty*nadj_;

  // Check if it is cheaper to calculate the full Jacobian and then multiply
  casadi_assert_warning(der_dir_cost <= 2*full_jac_cost, "Inefficient directional derivative calculation.");

  // Request an increase in the number of directional derivatives
  fcn_.requestNumSens(nfwd_,nadj_);
}

void DerivativeInternal::evaluate(int nfdir, int nadir){
  casadi_assert_message(nfdir==0, "Not implemeted");
  casadi_assert_message(nadir==0, "Not implemeted");
  
  casadi_log("DerivativeInternal::evaluate(" << nfdir << ", " << nadir<< "):begin  " << getOption("name"));
  
  // Number inputs and outputs of the function
  int num_in = fcn_.getNumInputs();
  int num_out = fcn_.getNumOutputs();
  
  // Number of derivative directions supported by the function
  int max_nfwd = fcn_->nfdir_;
  int max_nadj = fcn_->nadir_;
  
  // Current forward and adjoint direction
  int offset_fwd = 0, offset_adj = 0;
  
  // Current input/output index
  int iind=0;
  int oind=0;
  
  // Evaluate until everything has been determinated
  bool first_batch=true;
  while (first_batch || offset_fwd < nfwd_ || offset_adj < nadj_) {

    // Number of forward and adjoint directions in the current "batch"
    int nfwd_f_batch = std::min(nfwd_ - offset_fwd, max_nfwd);
    int nadj_f_batch = std::min(nadj_ - offset_adj, max_nadj);

    // Pass the argument to the function if first batch
    if(first_batch){
      for(int i=0; i<num_in; ++i){
        fcn_.setInput(input(iind++),i);
      }
    }
    
    // Pass the forward seeds to the function
    for(int d = 0; d < nfwd_f_batch; ++d){
      for(int i = 0; i < num_in; ++i){
        fcn_.setFwdSeed(input(iind++),i,d);
      }
    }

    // Pass the adjoint seed to the function
    for(int d = 0; d < nadj_f_batch; ++d){
      for(int i = 0; i < num_out; ++i) {
        fcn_.setAdjSeed(input(iind++),i,d);
      }
    }

    // Evaluate
    fcn_.evaluate(nfwd_f_batch, nadj_f_batch);
    
    // Get the outputs if first evaluation
    if(first_batch){
      for(int i = 0; i < num_out; ++i) {
        fcn_.getOutput(output(oind++),i);
      }
    }

    // Get the forward sensitivities
    for(int d = 0; d < nfwd_f_batch; ++d){
      for(int i = 0; i < num_out; ++i) {
        fcn_.getFwdSens(output(oind++),i,d);
      }
    }

    // Get the adjoint sensitivities
    for(int d = 0; d < nadj_f_batch; ++d){
      for(int i = 0; i < num_in; ++i) {
        fcn_.getAdjSens(output(oind++),i,d);
      }
    }

    // Update direction offsets
    offset_fwd += nfwd_f_batch;
    offset_adj += nadj_f_batch;
    
    // No longer first batch
    first_batch = false;
  }
  casadi_log("DerivativeInternal::evaluate(" << nfdir << ", " << nadir<< "):end  " << getOption("name"));
}

} // namespace CasADi

