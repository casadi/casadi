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

#include "external_function_internal.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

// The following works for Linux, something simular is needed for Windows
#ifdef WITH_JIT 
#include <dlfcn.h>
#endif // WITH_JIT 

namespace CasADi{

using namespace std;

ExternalFunctionInternal::ExternalFunctionInternal(const std::string& bin_name) : bin_name_(bin_name){
#ifdef WITH_JIT 

  // Load the dll
  handle_ = 0;
  handle_ = dlopen(bin_name_.c_str(), RTLD_LAZY);    
  if (!handle_) {
    stringstream ss;
    ss << "ExternalFunctionInternal: Cannot open function: " << bin_name_ << ". error code: "<< dlerror();
    throw CasadiException(ss.str());
  }
  dlerror(); // reset error

  // Initialize and get the number of inputs and outputs
  initPtr init = (initPtr)dlsym(handle_, "init");
  if(dlerror()) throw CasadiException("ExternalFunctionInternal: no \"init\" found");
  int n_in, n_out;
  int flag = init(&n_in, &n_out);
  if(flag) throw CasadiException("ExternalFunctionInternal: \"init\" failed");
  
  // Pass to casadi
  input_.resize(n_in);
  output_.resize(n_out);

  // Get the sparsity pattern
  getSparsityPtr getSparsity = (getSparsityPtr)dlsym(handle_, "getSparsity");
  if(dlerror()) throw CasadiException("ExternalFunctionInternal: no \"getSparsity\" found");

  for(int i=0; i<n_in+n_out; ++i){
    // Get sparsity from file
    int nrow, ncol, *rowind, *col;
    flag = getSparsity(i,&nrow,&ncol,&rowind,&col);
    if(flag) throw CasadiException("ExternalFunctionInternal: \"getSparsity\" failed");

    // Row offsets
    vector<int> rowindv(rowind,rowind+nrow+1);
    
    // Number of nonzeros
    int nnz = rowindv.back();
    
    // Columns
    vector<int> colv(col,col+nnz);
    
    // Sparsity
    CRSSparsity sp(nrow,ncol,colv,rowindv);
    
    // Save to inputs/outputs
    if(i<n_in){
      input(i) = Matrix<double>(sp,0);
    } else {
      output(i-n_in) = Matrix<double>(sp,0);
    }
  }
  
  //
  evaluate_ = (evaluatePtr) dlsym(handle_, "evaluate");
  if(dlerror()) throw CasadiException("ExternalFunctionInternal: no \"evaluate\" found");
  
#else // WITH_JIT 
  throw CasadiException("WITH_JIT  not activated");
#endif // WITH_JIT 
  
}
    
ExternalFunctionInternal* ExternalFunctionInternal::clone() const{
  throw CasadiException("Error ExternalFunctionInternal cannot be cloned");
}

ExternalFunctionInternal::~ExternalFunctionInternal(){
#ifdef WITH_JIT 
  // close the dll
  if(handle_) dlclose(handle_);
#endif // WITH_JIT 
}

void ExternalFunctionInternal::evaluate(int nfdir, int nadir){
#ifdef WITH_JIT 
  int flag = evaluate_(&input_array_[0],&output_array_[0]);
  if(flag) throw CasadiException("ExternalFunctionInternal: \"evaluate\" failed");
#endif // WITH_JIT 
}
  
void ExternalFunctionInternal::init(){
  // Call the init function of the base class
  FXInternal::init();

  // Get pointers to the inputs
  input_array_.resize(input_.size());
  for(int i=0; i<input_array_.size(); ++i)
    input_array_[i] = &input(i).front();

  // Get pointers to the outputs
  output_array_.resize(output_.size());
  for(int i=0; i<output_array_.size(); ++i)
    output_array_[i] = &output(i).front();
}




} // namespace CasADi

