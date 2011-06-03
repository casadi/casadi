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

#include "parallelizer_internal.hpp"
#include "mx_function.hpp"

using namespace std;

namespace CasADi{
  
ParallelizerInternal::ParallelizerInternal(const std::vector<FX>& funcs) : funcs_(funcs){
  addOption("parallelization", OT_STRING, "serial"); // serial, openmp or mpi
  addOption("save_corrected_input", OT_BOOLEAN, false);
}

ParallelizerInternal::~ParallelizerInternal(){
}


void ParallelizerInternal::init(){
  // Get mode
  if(getOption("parallelization")=="serial")
    mode_ = SERIAL;
  else if(getOption("parallelization")=="openmp")
    mode_ = OPENMP;
  else if(getOption("parallelization")=="mpi")
    mode_ = MPI;
  else
    throw CasadiException(string("Parallelization mode: ")+getOption("parallelization").toString());

  // Check if a node is a copy of another
  copy_of_.resize(funcs_.size(),-1);
  map<void*,int> is_copy_of;
  for(int i=0; i<funcs_.size(); ++i){
    // Check if the function has already been assigned an index
    map<void*,int>::const_iterator it=is_copy_of.find(funcs_[i].get());
    if(it!=is_copy_of.end()){
      copy_of_[i] = it->second;
    } else {
      is_copy_of[funcs_[i].get()] = i;
    }
  }
  
  // Initialize the dependend functions
  for(vector<FX>::iterator it=funcs_.begin(); it!=funcs_.end(); ++it){
    // Make sure that the functions are unique if we are using OpenMP
    if(mode_==OPENMP)
      it->makeUnique();
    
    // Initialize
    it->init();
  }
  
  // Clear the indices
  inind_.clear();   inind_.push_back(0);
  outind_.clear();  outind_.push_back(0);
  
  // Add the inputs and outputs
  for(vector<FX>::iterator it=funcs_.begin(); it!=funcs_.end(); ++it){
    inind_.push_back(inind_.back()+it->getNumInputs());
    outind_.push_back(outind_.back()+it->getNumOutputs());
  }
  setNumInputs(inind_.back());
  setNumOutputs(outind_.back());
  
  // Copy inputs and output dimensions and structure
  for(int i=0; i<funcs_.size(); ++i){
    for(int j=inind_[i]; j<inind_[i+1]; ++j)
      input(j) = funcs_[i].input(j-inind_[i]);
    for(int j=outind_[i]; j<outind_[i+1]; ++j)
      output(j) = funcs_[i].output(j-outind_[i]);
  }
  
  // Call the init function of the base class
  FXInternal::init();
  
  // Should corrected input values be saved after evaluation?
  save_corrected_input_ = getOption("save_corrected_input").toInt();
}

void ParallelizerInternal::evaluate(int nfdir, int nadir){
  switch(mode_){
    case SERIAL:
    {
      for(int task=0; task<funcs_.size(); ++task)
        evaluateTask(task,nfdir,nadir);
      break;
    }
    case OPENMP:
    {
      #pragma omp parallel for
      for(int task=0; task<funcs_.size(); ++task)
        evaluateTask(task,nfdir,nadir);
      break;
    }
    case MPI:
    {
      throw CasadiException("ParallelizerInternal::evaluate: MPI not implemented");
      break;
    }
  }
}

void ParallelizerInternal::evaluateTask(int task, int nfdir, int nadir){
  
  // Get a reference to the function
  FX& fcn = funcs_[task];
  
  // Copy inputs to functions
  for(int j=inind_[task]; j<inind_[task+1]; ++j){
    fcn.input(j-inind_[task]).set(input(j));
  }
  
  // Copy forward seeds
  for(int dir=0; dir<nfdir; ++dir){
    for(int j=inind_[task]; j<inind_[task+1]; ++j){
      fcn.fwdSeed(j-inind_[task],dir).set(fwdSeed(j,dir));
    }
  }
  
  // Copy adjoint seeds
  for(int dir=0; dir<nadir; ++dir){
    for(int j=outind_[task]; j<outind_[task+1]; ++j){
      fcn.adjSeed(j-outind_[task],dir).set(adjSeed(j,dir));
    }
  }

  // Evaluate
  fcn.evaluate(nfdir, nadir);
    
  // Get the results
  for(int j=outind_[task]; j<outind_[task+1]; ++j){
    fcn.output(j-outind_[task]).get(output(j));
  }
  
  // Get the forward sensitivities
  for(int dir=0; dir<nfdir; ++dir){
    for(int j=outind_[task]; j<outind_[task+1]; ++j){
      fcn.fwdSens(j-outind_[task],dir).get(fwdSens(j,dir));
    }
  }

  // Get the adjoint sensitivities
  for(int dir=0; dir<nadir; ++dir){
    for(int j=inind_[task]; j<inind_[task+1]; ++j){
      fcn.adjSens(j-inind_[task],dir).get(adjSens(j,dir));
    }
  }
  
  // Save corrected input values // TODO: REMOVE!
  if(save_corrected_input_){
    for(int j=inind_[task]; j<inind_[task+1]; ++j){
      fcn.getInput(input(j),j-inind_[task]);
    }
  }
}

FX ParallelizerInternal::jacobian(int iind, int oind){
  // Number of tasks
  int ntask = inind_.size()-1;
  
  // Find out which task corresponds to the iind
  int task;
  for(task=0; task<ntask && iind>=inind_[task+1]; ++task);
  
  // Check if the output index is also in this task
  if(oind>=outind_[task] && oind<outind_[task+1]){
    // Call the base class
/*    return funcs_[task].jacobian(iind-inind_[task],oind-outind_[task]);*/
    return FXInternal::jacobian(iind,oind);
  } else {
    // Return a null reference
    return FX();
  }
}

CRSSparsity ParallelizerInternal::getJacSparsity(int iind, int oind){
  // Number of tasks
  int ntask = inind_.size()-1;
  
  // Find out which task corresponds to the iind
  int task;
  for(task=0; task<ntask && iind>=inind_[task+1]; ++task);
  
  // Check if the output index is also in this task
  if(oind>=outind_[task] && oind<outind_[task+1]){
    // Call the base class
    return FXInternal::getJacSparsity(iind,oind);
  } else {
    // All-zero jacobian
    return CRSSparsity();
  }
}

FX ParallelizerInternal::jacobian(const vector<pair<int,int> >& jblocks){
  
  // Symbolic input
  vector<MX> j_in(getNumInputs());
  for(int i=0; i<j_in.size(); ++i){
    stringstream name;
    name << "x_" << i;
    j_in[i] = MX(name.str(),input(i).sparsity());
  }
  
  // Jacobian functions to be evaluated in parallel
  vector<FX> jac_funcs(funcs_.size());
  
  // Current jacobian block
  vector<pair<int,int> >::const_iterator jit = jblocks.begin();
  
  // Loop over tasks
  for(int i=0; i<funcs_.size(); ++i){
    // Local jacobian blocks
    vector<pair<int,int> > jblocks_local;

    // Loop over jacobian blocks
    while(jit->first >= outind_[i] && jit->first < outind_[i+1] && jit->second < inind_[i+1] && (jit->second == -1 || jit->second >= inind_[i])){
      jblocks_local.push_back(pair<int,int>(jit->first - outind_[i], jit->second == -1 ? -1 : jit->second - inind_[i]));
      jit++;
    }
    
    // Check if the function has already been calculated
    if(copy_of_[i]>=0){
      // The Jacobian has already been calculated
      jac_funcs[i] = jac_funcs[copy_of_[i]];
    } else {
      // Differentiate the function
      jac_funcs[i] = funcs_[i].jacobian(jblocks_local);
    }
  }

  // Make sure that all the blocks have been visited
  casadi_assert(jit == jblocks.end());

  // Create new Parallelizer for the functions and Jacobians
  Parallelizer parjac(jac_funcs);
  parjac.setOption(dictionary()); // copy the options
  
  return parjac;
}


} // namespace CasADi

