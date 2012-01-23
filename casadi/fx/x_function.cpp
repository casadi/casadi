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

vector<SXMatrix> XFunction::evalSX(const vector<SXMatrix>& arg){
  casadi_assert_message(isInit(),"Function has not been initialized");
  
  // Copy the arguments into a new vector with the right sparsity
  casadi_assert(arg.size()==getNumInputs());
  vector<SXMatrix> arg2 = arg;
  for(int iind=0; iind<arg.size(); ++iind){
    // If sparsities do not match, we need to map the nonzeros
    if(!(arg2[iind].sparsity()==input(iind).sparsity())){
      // The sparsity should be that of the inputs
      arg2[iind] = SXMatrix(input(iind).sparsity(),0);

      // Make sure that the dimensions match
      casadi_assert(arg[iind].size1()==arg2[iind].size1());
      casadi_assert(arg[iind].size2()==arg2[iind].size2());

      // Get the indices of the known supplied arguments
      vector<int> known_ind = arg[iind].sparsity().getElementMapping(false);
      
      // Find the corresponding nonzeros of the argument matrix
      arg2[iind].sparsity().getNZInplace(known_ind);
      
      // Set the element values
      for(int k=0; k<known_ind.size(); ++k){
        if(known_ind[k]!=-1){
          arg2[iind].at(known_ind[k]) = arg[iind].at(k);
        }
      }
    }
  }

  // Create result vector with correct sparsity for the result
  vector<SXMatrix> res(getNumOutputs());
  for(int i=0; i<res.size(); ++i){
    res[i] = SXMatrix(output(i).sparsity());
  }
  
  // No sensitivities
  vector<vector<SXMatrix> > dummy;
  
  // Evaluate the algorithm
  (*this)->eval(arg2,res,dummy,dummy,dummy,dummy,false,false);
  
  // Return the result
  return res;
}

vector<MX> XFunction::evalMX(const vector<MX>& arg){
  vector<MX> res;
  vector<vector<MX> > dummy;
  (*this)->evalMX(arg,res,dummy,dummy,dummy,dummy,false,false);
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

