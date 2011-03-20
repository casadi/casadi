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

void ParallelizerInternal::evaluate_new(int nfdir, int nadir){
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
  fcn.evaluate_new(nfdir, nadir);
    
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


} // namespace CasADi

