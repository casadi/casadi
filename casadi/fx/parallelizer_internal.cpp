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
  addOption("mode", OT_STRING, "serial"); // serial, openmp or mpi
}

ParallelizerInternal::~ParallelizerInternal(){
}


void ParallelizerInternal::init(){
  // Initialize the dependend functions
  for(vector<FX>::iterator it=funcs_.begin(); it!=funcs_.end(); ++it)
    it->init();
  
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

  // Get mode
  if(getOption("mode")=="serial")
    mode_ = SERIAL;
  else if(getOption("mode")=="openmp")
    mode_ = OPENMP;
  else if(getOption("mode")=="mpi")
    mode_ = MPI;
  else
    throw CasadiException(string("Unknown mode: ")+getOption("mode").toString());
}

void ParallelizerInternal::evaluate(int fsens_order, int asens_order){
  switch(mode_){
    case SERIAL:
    {
      for(int task=0; task<funcs_.size(); ++task)
        evaluateTask(task,fsens_order,asens_order);
      break;
    }
    case OPENMP:
    {
      #pragma omp parallel for
      for(int task=0; task<funcs_.size(); ++task)
        evaluateTask(task,fsens_order,asens_order);
      break;
    }
    case MPI:
    {
      throw CasadiException("ParallelizerInternal::evaluate: MPI not implemented");
      break;
    }
  }
}

void ParallelizerInternal::evaluateTask(int task, int fsens_order, int asens_order){
  // Get a reference to the function
  FX& fcn = funcs_[task];
  
  // Copy inputs to functions TODO: add seeds
  for(int j=inind_[task]; j<inind_[task+1]; ++j)
    fcn.input(j-inind_[task]).set(input(j));

  // Evaluate
  fcn.evaluate(fsens_order,asens_order);
    
  // Get the results TODO: add derivatives
  for(int j=outind_[task]; j<outind_[task+1]; ++j)
    fcn.output(j-outind_[task]).get(output(j));
}


} // namespace CasADi

