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

#include "mx_function_internal.hpp"
#include "../mx/mx_node.hpp"
#include "../stl_vector_tools.hpp"

#include <stack>
#include <typeinfo>
#include <cassert>

using namespace std;

namespace CasADi{

bool MXFunction::checkNode() const{
  return dynamic_cast<const MXFunctionInternal*>(get())!=0;
}

MXFunction::MXFunction(){
}

MXFunction::MXFunction(const MX& inputm, const MX& outputm){
  vector<MX> inputv(1);
  inputv[0] = inputm;
  vector<MX> outputv(1);
  outputv[0] = outputm;
  assignNode(new MXFunctionInternal(inputv,outputv));
}

MXFunction::MXFunction(const MX& inputm, const std::vector<MX>& outputv){
  vector<MX> inputv(1);
  inputv[0] = inputm;
  assignNode(new MXFunctionInternal(inputv,outputv));
}

MXFunction::MXFunction(const std::vector<MX>& inputv, const MX& outputm){
  vector<MX> outputv(1);
  outputv[0] = outputm;
  assignNode(new MXFunctionInternal(inputv,outputv));
}

MXFunction::MXFunction(const std::vector<MX>& inputv, const std::vector<MX>& outputv){
  assignNode(new MXFunctionInternal(inputv,outputv));
}

const MXFunctionInternal* MXFunction::operator->() const{
  return (const MXFunctionInternal*)FX::operator->();
}

MXFunctionInternal* MXFunction::operator->(){
  return (MXFunctionInternal*)FX::operator->();
}

const MX MXFunction::inputMX(int iind) const{
  return (*this)->inputv.at(iind);
}

const MX MXFunction::outputMX(int oind) const{
  return (*this)->outputv.at(oind);
}

const std::vector<MXAlgEl>& MXFunction::algorithm() const{
  return (*this)->alg;
}

int MXFunction::countNodes() const{
  casadi_assert(isInit());
  return algorithm().size();
}

void MXFunction::setLiftingFunction(LiftingFunction liftfun, void* user_data){
  (*this)->setLiftingFunction(liftfun,user_data);
}

std::vector<MX> MXFunction::jac(int iind){
  return (*this)->jac(iind);
}

std::vector<MX> MXFunction::grad(int oind){
  return (*this)->grad(oind);
}

SXFunction MXFunction::expand(const std::vector<SXMatrix>& inputv){
  return (*this)->expand(inputv);
}


} // namespace CasADi

