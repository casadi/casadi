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

#include "functor_internal.hpp"

#include "matrix/crs_sparsity.hpp"
#include "fx/fx.hpp"
#include "fx/custom_function.hpp"

using namespace std;

namespace CasADi{

  //Call

  SparsityGeneratorCInternal::SparsityGeneratorCInternal(SparsityGeneratorCPtr ptr)  : FunctorCInternal(ptr) {
  }
  
  CRSSparsity SparsityGeneratorCInternal::call(FX& fcn, int iind, int oind, void* user_data) {
    casadi_assert(ptr_!=0);
    return ptr_(fcn, iind, oind, user_data);
  }
  
  SparsityGeneratorCInternal* SparsityGeneratorCInternal::clone() const {
    return new SparsityGeneratorCInternal(ptr_);
  }
  
  JacobianGeneratorCInternal::JacobianGeneratorCInternal(JacobianGeneratorCPtr ptr)  : FunctorCInternal(ptr) {
  }
  
  FX JacobianGeneratorCInternal::call(FX& fcn, int iind, int oind, void* user_data) {
    casadi_assert(ptr_!=0);
    return ptr_(fcn, iind, oind, user_data);
  }
  
  JacobianGeneratorCInternal* JacobianGeneratorCInternal::clone() const {
    return new JacobianGeneratorCInternal(ptr_);
  }
  
  CustomEvaluateCInternal::CustomEvaluateCInternal(CustomEvaluateCPtr ptr)  : FunctorCInternal(ptr) {
  }
  
  void CustomEvaluateCInternal::call(CustomFunction& fcn, int iind, int oind, void* user_data) {
    casadi_assert(ptr_!=0);
    ptr_(fcn, iind, oind, user_data);
  }
  
  CustomEvaluateCInternal* CustomEvaluateCInternal::clone() const {
    return new CustomEvaluateCInternal(ptr_);
  }

  CallbackCInternal::CallbackCInternal(CallbackCPtr ptr)  : FunctorCInternal(ptr) {
  }
  
  int CallbackCInternal::call(FX& fcn, void* user_data) {
    casadi_assert(ptr_!=0);
    return ptr_(fcn, user_data);
  }
  
  CallbackCInternal* CallbackCInternal::clone() const {
    return new CallbackCInternal(ptr_);
  }
  
} // namespace CasADi

