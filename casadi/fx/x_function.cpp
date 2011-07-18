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

#include "x_function_internal.hpp"

namespace CasADi{

using namespace std;


XFunction::XFunction(){
}

const XFunctionInternal* XFunction::operator->() const{
  return (const XFunctionInternal*)FX::operator->();
}

XFunctionInternal* XFunction::operator->(){
  return (XFunctionInternal*)FX::operator->();
}

bool XFunction::checkNode() const{
  return dynamic_cast<const XFunctionInternal*>(get())!=0;
}

vector<SXMatrix> XFunction::eval(const vector<SXMatrix>& arg){
  casadi_assert_message(isInit(),"Function has not been initialized");

  // Create result vector with correct sparsity for the result
  vector<SXMatrix> res(getNumOutputs());
  for(int i=0; i<res.size(); ++i){
    res[i] = SXMatrix(output(i).sparsity());
  }
  
  // Evaluate the algorithm
  (*this)->evaluateSX(arg,res);
  
  // Return the result
  return res;
}

SXMatrix XFunction::eval(const SXMatrix& arg){
  return eval(vector<SXMatrix>(1,arg)).at(0);
}

vector< vector<SX> > XFunction::eval(const vector< vector<SX> >& arg){
  // Convert input
  casadi_assert(getNumInputs()==arg.size());
  vector<SXMatrix> argv(getNumInputs());
  for(int i=0; i<arg.size(); ++i){
    casadi_assert(arg[i].size()==input(i).size());
    argv[i] = SXMatrix(input(i).sparsity(),arg[i]);
  }
  
  // Call the normal eval function
  vector<SXMatrix> resv = eval(argv);
  
  // Convert result to a vector of vectors
  vector< vector<SX> > res(resv.size());
  for(int i=0; i<res.size(); ++i){
    res[i] = resv[i].data();
  }
  
  return res;
}

vector<SX> XFunction::eval(const vector<SX>& arg){
  return eval(vector< vector<SX> >(1,arg)).at(0);
}


} // namespace CasADi

