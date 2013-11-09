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
  

CRSSparsity SparsityGenerator::operator()(FX& fcn, int iind, int oind, void* user_data) {
  checkNode();
  return static_cast<SparsityGeneratorInternal*>(SharedObject::operator->())->call(fcn, iind, oind, user_data);
}

FX JacobianGenerator::operator()(FX& fcn, int iind, int oind, void* user_data) {
  checkNode();
  return static_cast<JacobianGeneratorInternal*>(SharedObject::operator->())->call(fcn, iind, oind, user_data);
}
  
void CustomEvaluate::operator()(CustomFunction& fcn, int nfdir, int nadir, void* user_data) {
  checkNode();
  static_cast<CustomEvaluateInternal*>(SharedObject::operator->())->call(fcn, nfdir, nadir, user_data);
}
  
int Callback::operator()(FX& fcn, void* user_data) {
  checkNode();
  return static_cast<CallbackInternal*>(SharedObject::operator->())->call(fcn, user_data);
}
  
SparsityGeneratorC::SparsityGeneratorC(SparsityGeneratorCPtr ptr) {
  assignNode(new SparsityGeneratorCInternal(ptr));
}
  
JacobianGeneratorC::JacobianGeneratorC(JacobianGeneratorCPtr ptr) {
  assignNode(new JacobianGeneratorCInternal(ptr));
}  

CustomEvaluateC::CustomEvaluateC(CustomEvaluateCPtr ptr) {
  assignNode(new CustomEvaluateCInternal(ptr));
}

CallbackC::CallbackC(CallbackCPtr ptr) {
  assignNode(new CallbackCInternal(ptr));
}
  
} // namespace CasADi
