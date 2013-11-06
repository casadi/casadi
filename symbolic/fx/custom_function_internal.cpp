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

#include "custom_function_internal.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

namespace CasADi{

using namespace std;

CustomFunctionInternal::CustomFunctionInternal(const CustomEvaluate &c_fcn, const std::vector<CasADi::CRSSparsity> &inputscheme, const std::vector<CasADi::CRSSparsity> &outputscheme) : evaluate_(c_fcn){
  setNumInputs(inputscheme.size());
  setNumOutputs(outputscheme.size());
  
  for (int k=0;k<inputscheme.size();k++) {
    input(k) = DMatrix(inputscheme[k],0);
  }
  
  for (int k=0;k<outputscheme.size();k++) {
    output(k) = DMatrix(outputscheme[k],0);
  }
  
  // Make the ref object a non-refence counted pointer to this (as reference counting would prevent deletion of the object)
  ref_.assignNodeNoCount(this);
  
  setOption("max_number_of_fwd_dir",0);
  setOption("max_number_of_adj_dir",0);
  
}

CustomFunctionInternal::~CustomFunctionInternal(){
  // Explicitly remove the pointer to this (as the counter would otherwise be decreased)
  ref_.assignNodeNoCount(0);
}

void CustomFunctionInternal::evaluate(int nfdir, int nadir){
  casadi_assert_message(!evaluate_.isNull(), "CustomFunctionInternal::evaluate: pointer is null");
  evaluate_(ref_,nfdir,nadir,user_data_);  
}

void CustomFunctionInternal::init(){
  FXInternal::init();
}


} // namespace CasADi

