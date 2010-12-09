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

#include "external_function.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <dlfcn.h>
#include <cassert>

namespace CasADi{

using namespace std;

#if 0  

ExternalFunctionNode::ExternalFunctionNode(const string& bin_name){
  int flag;

  // Load the dll
  handle = dlopen(bin_name.c_str(), RTLD_LAZY);    
  if (!handle) {
        cerr << "Cannot open function: " << bin_name << ". error code: "<< dlerror() << endl;
        throw "ExternalFunctionNode::ExternalFunctionNode";
  }
  dlerror(); // reset error

  // Initialize and get the number of inputs and outputs
  typedef int (*initPtr)(int*,int*);
  initPtr init = (initPtr)dlsym(handle, "init");
  assert(!dlerror()); // make sure that the function was found
  int n_in, n_out;
  flag = init(&n_in, &n_out);
  assert(flag==0);
  
  // Get the size of the inputs
  typedef int (*getInputSizePtr)(int,int*,int*);
  getInputSizePtr getInputSize = (getInputSizePtr)dlsym(handle, "getInputSize");
  assert(!dlerror()); // make sure that the function was found
  
  input.resize(n_in);
  for(int i=0; i<n_in; ++i){
    MatrixSize sz;
    flag = getInputSize(i,&sz.nrow,&sz.ncol);
    input[i] = FunctionIO(sz);
    assert(flag==0);
  }

  // Get the size of the outputs
  typedef int (*getOutputSizePtr)(int,int*,int*);
  getOutputSizePtr getOutputSize = (getOutputSizePtr)dlsym(handle, "getOutputSize");
  assert(!dlerror()); // make sure that the function was found
  output.resize(n_out);
  for(int i=0; i<n_out; ++i){
    MatrixSize sz;
    flag = getOutputSize(i,&sz.nrow,&sz.ncol);
    output[i] = FunctionIO(sz);
    assert(flag==0);
  }
    
  // Try to read the functions and set to null if not found
  setArgument_ptr = (setterPtr) dlsym(handle, "setArgument");
  if(dlerror()) setArgument_ptr = 0;

  getArgument_ptr = (getterPtr) dlsym(handle, "getArgument");
  if(dlerror()) getArgument_ptr = 0;

  setResult_ptr = (setterPtr) dlsym(handle, "setResult");
  if(dlerror()) setResult_ptr = 0;

  getResult_ptr = (getterPtr) dlsym(handle, "getResult");
  if(dlerror()) getResult_ptr = 0;

  clear_ptr = (clearerPtr) dlsym(handle, "clear");
  if(dlerror()) clear_ptr = 0;

  evaluate_ptr = (evaluaterPtr) dlsym(handle, "evaluate");
  if(dlerror()) evaluate_ptr = 0;

  evaluateFwd_ptr = (evaluaterPtr) dlsym(handle, "evaluateFwd");
  if(dlerror()) evaluateFwd_ptr = 0;

  evaluateAdj_ptr = (evaluaterPtr) dlsym(handle, "evaluateAdj");
  if(dlerror()) evaluateAdj_ptr = 0;
  

}

ExternalFunctionNode::~ExternalFunctionNode(){
    int flag;

    // close the dll
    dlclose(handle);
    
  }

void ExternalFunctionNode::setArgument(const double *x, int ind, int ord){
  assert(setArgument_ptr);
  int flag = setArgument_ptr(x,ind,ord); 
  assert(flag==0);
}

void ExternalFunctionNode::getArgument(double *x, int ind, int ord) const{
  assert(getArgument_ptr);
  int flag = getArgument_ptr(x,ind,ord); 
  assert(flag==0);
}

void ExternalFunctionNode::setResult(const double *x, int ind, int ord){ 
  assert(setResult_ptr);
  int flag = setResult_ptr(x,ind,ord);
  assert(flag == 0);
}


void ExternalFunctionNode::getResult(double *x, int ind, int ord) const{
  assert(getResult_ptr);
  int flag = getResult_ptr(x,ind,ord);
  assert(flag == 0);
}

void ExternalFunctionNode::clear(int ord){ 
  assert(clear_ptr);
  int flag = clear_ptr(ord);
  assert(flag == 0);
}

void ExternalFunctionNode::evaluate(int tape_order){ 
  assert(evaluate_ptr);
  int flag = evaluate_ptr();
  assert(flag == 0);
}

void ExternalFunctionNode::evaluateFwd(bool use_tape){ 
  assert(evaluateFwd_ptr);
  int flag = evaluateFwd_ptr();
  assert(flag == 0);
}

void ExternalFunctionNode::evaluateAdj(){ 
  assert(evaluateAdj_ptr);
  int flag = evaluateAdj_ptr();
  assert(flag == 0);
}

ExternalFunction::ExternalFunction(const string& bin_name){
  (FX&)(*this) = FX(new ExternalFunctionNode(bin_name));
}

ExternalFunctionNode* const ExternalFunction::get() const{
  return (ExternalFunctionNode*)node;
}

const ExternalFunctionNode* ExternalFunction::operator->() const{
  return (const ExternalFunctionNode*)node;
}

ExternalFunctionNode* ExternalFunction::operator->(){
  return (ExternalFunctionNode*)node;
}

ExternalFunction& ExternalFunction::operator=(const ExternalFunction& fx){
  (FX&)(*this) = fx;
}

void ExternalFunctionNode::init(){
  // Call the init function of the base class
  FXNode::init();
}

#endif


} // namespace CasADi

