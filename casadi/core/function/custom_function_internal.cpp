/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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

namespace casadi {

using namespace std;

CustomFunctionInternal::CustomFunctionInternal(
    const CustomEvaluate &c_fcn,
    const std::vector<casadi::Sparsity> &inputscheme,
    const std::vector<casadi::Sparsity> &outputscheme) : evaluate_(c_fcn) {
  setNumInputs(inputscheme.size());
  setNumOutputs(outputscheme.size());

  for (int k=0;k<inputscheme.size();k++) {
    input(k) = DMatrix(inputscheme[k], 0);
  }

  for (int k=0;k<outputscheme.size();k++) {
    output(k) = DMatrix(outputscheme[k], 0);
  }

  // Make the ref object a non-refence counted pointer to this (as reference counting
  // would prevent deletion of the object)
  ref_.assignNodeNoCount(this);

}

CustomFunctionInternal::~CustomFunctionInternal() {
  // Explicitly remove the pointer to this (as the counter would otherwise be decreased)
  ref_.assignNodeNoCount(0);
}

void CustomFunctionInternal::evaluate() {
  casadi_assert_message(!evaluate_.isNull(), "CustomFunctionInternal::evaluate: pointer is null");
  evaluate_(ref_, user_data_);
}

void CustomFunctionInternal::init() {
  FunctionInternal::init();
}


} // namespace casadi

