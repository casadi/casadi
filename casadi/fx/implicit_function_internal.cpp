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

#include "implicit_function_internal.hpp"

using namespace std;
namespace CasADi{

ImplicitFunctionInternal::ImplicitFunctionInternal(const FX& f, int nrhs) : f_(f), nrhs_(nrhs){
}

void ImplicitFunctionInternal::init(){
  // Initialize the residual function
  f_.init();
  
  // Allocate inputs
  setNumInputs(f_.getNumInputs()-1);
  for(int i=0; i<getNumInputs(); ++i)
    input(i) = f_.input(i+1);
  
  // Allocate outputs
  setNumOutputs(1);
  output() = f_.input(0);

  // Call the base class initializer
  FXInternal::init();

  // Number of equations
  N_ = output().size();

  // Get the number of directions of the function
  nfdir_fcn_ = f_.getOption("number_of_fwd_dir");
  nadir_fcn_ = f_.getOption("number_of_adj_dir");

}

ImplicitFunctionInternal::~ImplicitFunctionInternal(){
}
 
 
} // namespace CasADi

  


