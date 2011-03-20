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

#include "c_function_internal.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

namespace CasADi{

using namespace std;

CFunctionInternal::CFunctionInternal(CFunctionWrapper c_fcn) : evaluate_(c_fcn){
  user_data_ = 0;
  
  // Make the ref object a non-refence counted pointer to this (as reference counting would prevent deletion of the object)
  ref_.assignNodeNoCount(this);
}

CFunctionInternal::~CFunctionInternal(){
  // Explicitly remove the pointer to this (as the counter would otherwise be decreased)
  ref_.assignNodeNoCount(0);
}

void CFunctionInternal::setUserData(void* user_data){
  user_data_ = user_data;
}

void CFunctionInternal::evaluate_new(int nfdir, int nadir){
  casadi_assert_message(evaluate_!=0, "CFunctionInternal::evaluate: pointer is null");
  evaluate_(ref_,nfdir,nadir,user_data_);  
}

void CFunctionInternal::init(){
  FXInternal::init();
}


} // namespace CasADi

