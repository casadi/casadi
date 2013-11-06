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

#include "compiled_function_internal.hpp"
#include "../stl_vector_tools.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

namespace CasADi{

using namespace std;

CompiledFunctionInternal::CompiledFunctionInternal(const std::string& bin_name) : bin_name_(bin_name){
#ifdef WITH_DL 

  // Load the dll
#ifdef _WIN32
  handle_ = LoadLibrary(TEXT(bin_name_.c_str()));  
  casadi_assert_message(handle_!=0,"CompiledFunctionInternal: Cannot open function: " << bin_name_ << ". error code (WIN32): "<< GetLastError());

  initPtr init = (initPtr)GetProcAddress(handle_,TEXT("init"));
  if(init==0) throw CasadiException("CompiledFunctionInternal: no \"init\" found");
  getSparsityPtr getSparsity = (getSparsityPtr)GetProcAddress(handle_, TEXT("getSparsity"));
  if(getSparsity==0) throw CasadiException("CompiledFunctionInternal: no \"getSparsity\" found");
  evaluate_ = (evaluatePtr) GetProcAddress(handle_, TEXT("evaluateWrap"));
  if(evaluate_==0) throw CasadiException("CompiledFunctionInternal: no \"evaluateWrap\" found");

#else // _WIN32
  handle_ = dlopen(bin_name_.c_str(), RTLD_LAZY);  
  casadi_assert_message(handle_!=0,"CompiledFunctionInternal: Cannot open function: " << bin_name_ << ". error code: "<< dlerror());

  // reset error
  dlerror(); 

  // Load symbols
  initPtr init = (initPtr)dlsym(handle_, "init");
  if(dlerror()) throw CasadiException("CompiledFunctionInternal: no \"init\" found");
  getSparsityPtr getSparsity = (getSparsityPtr)dlsym(handle_, "getSparsity");
  if(dlerror()) throw CasadiException("CompiledFunctionInternal: no \"getSparsity\" found");
  evaluate_ = (evaluatePtr) dlsym(handle_, "evaluateWrap");
  if(dlerror()) throw CasadiException("CompiledFunctionInternal: no \"evaluateWrap\" found");
#endif // _WIN32

  // Initialize and get the number of inputs and outputs
  int n_in=-1, n_out=-1;
  int flag = init(&n_in, &n_out);
  if(flag) throw CasadiException("CompiledFunctionInternal: \"init\" failed");
  
  // Pass to casadi
  input_.resize(n_in);
  output_.resize(n_out);
  
  // Get the sparsity pattern
  for(int i=0; i<n_in+n_out; ++i){
    // Get sparsity from file
    int nrow, ncol, *rowind, *col;
    flag = getSparsity(i,&nrow,&ncol,&rowind,&col);
    if(flag) throw CasadiException("CompiledFunctionInternal: \"getSparsity\" failed");

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
    
#else // WITH_DL 
  throw CasadiException("WITH_DL  not activated");
#endif // WITH_DL 
  
}
    
CompiledFunctionInternal* CompiledFunctionInternal::clone() const{
  throw CasadiException("Error CompiledFunctionInternal cannot be cloned");
}

CompiledFunctionInternal::~CompiledFunctionInternal(){
#ifdef WITH_DL 
  // close the dll
#ifdef _WIN32
  if(handle_) FreeLibrary(handle_);
#else // _WIN32
  if(handle_) dlclose(handle_);
#endif // _WIN32
#endif // WITH_DL 
}

void CompiledFunctionInternal::evaluate(int nfdir, int nadir){
#ifdef WITH_DL 
  int flag = evaluate_(getPtr(input_array_),getPtr(output_array_));
  if(flag) throw CasadiException("CompiledFunctionInternal: \"evaluate\" failed");
#endif // WITH_DL 
}
  
void CompiledFunctionInternal::init(){
  // Call the init function of the base class
  FXInternal::init();

  // Get pointers to the inputs
  input_array_.resize(input_.size());
  for(int i=0; i<input_array_.size(); ++i)
    input_array_[i] = input(i).ptr();

  // Get pointers to the outputs
  output_array_.resize(output_.size());
  for(int i=0; i<output_array_.size(); ++i)
    output_array_[i] = output(i).ptr();
}




} // namespace CasADi

