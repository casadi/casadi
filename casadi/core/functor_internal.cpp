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


#include "functor_internal.hpp"

#include "matrix/sparsity.hpp"
#include "function/function.hpp"
#include "function/custom_function.hpp"

using namespace std;

namespace casadi {

  //Call

  DerivativeGeneratorCInternal::DerivativeGeneratorCInternal(DerivativeGeneratorCPtr ptr) :
      FunctorCInternal<DerivativeGeneratorCPtr>(ptr) {
  }

  Function DerivativeGeneratorCInternal::call(Function& fcn, int nfwd, int nadj, void* user_data) {
    casadi_assert(ptr_!=0);
    return ptr_(fcn, nfwd, nadj, user_data);
  }

  DerivativeGeneratorCInternal* DerivativeGeneratorCInternal::clone() const {
    return new DerivativeGeneratorCInternal(ptr_);
  }

  CustomEvaluateCInternal::CustomEvaluateCInternal(CustomEvaluateCPtr ptr) :
      FunctorCInternal<CustomEvaluateCPtr>(ptr) {
  }

  void CustomEvaluateCInternal::call(CustomFunction& fcn, void* user_data) {
    casadi_assert(ptr_!=0);
    ptr_(fcn, user_data);
  }

  CustomEvaluateCInternal* CustomEvaluateCInternal::clone() const {
    return new CustomEvaluateCInternal(ptr_);
  }

  CallbackCInternal::CallbackCInternal(CallbackCPtr ptr)  : FunctorCInternal<CallbackCPtr>(ptr) {
  }

  int CallbackCInternal::call(Function& fcn, void* user_data) {
    casadi_assert(ptr_!=0);
    return ptr_(fcn, user_data);
  }

  CallbackCInternal* CallbackCInternal::clone() const {
    return new CallbackCInternal(ptr_);
  }

} // namespace casadi
