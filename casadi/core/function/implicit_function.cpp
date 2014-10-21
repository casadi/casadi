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


#include "implicit_function_internal.hpp"

using namespace std;
namespace casadi {

  ImplicitFunction::ImplicitFunction() {
  }

  ImplicitFunction::ImplicitFunction(const std::string& name, const Function& f,
                                     const Function& jac,
                                     const LinearSolver& linsol) {
    assignNode(ImplicitFunctionInternal::instantiatePlugin(name, f));
    setOption("jacobian_function", jac);
    setOption("linear_solver_function", linsol);
  }

  ImplicitFunctionInternal* ImplicitFunction::operator->() {
    return static_cast<ImplicitFunctionInternal*>(Function::operator->());
  }

  const ImplicitFunctionInternal* ImplicitFunction::operator->() const {
    return static_cast<const ImplicitFunctionInternal*>(Function::operator->());
  }

  bool ImplicitFunction::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const ImplicitFunctionInternal*>(ptr)!=0;
  }

  Function& ImplicitFunction::getF() {
    casadi_assert(!isNull());
    return (*this)->f_;
  }

  Function& ImplicitFunction::getJac() {
    casadi_assert(!isNull());
    return (*this)->jac_;
  }

  LinearSolver& ImplicitFunction::getLinsol() {
    casadi_assert(!isNull());
    return (*this)->linsol_;
  }

  bool ImplicitFunction::hasPlugin(const std::string& name) {
    return ImplicitFunctionInternal::hasPlugin(name);
  }

  void ImplicitFunction::loadPlugin(const std::string& name) {
    ImplicitFunctionInternal::loadPlugin(name);
  }

  std::string ImplicitFunction::doc(const std::string& name) {
    return ImplicitFunctionInternal::getPlugin(name).doc;
  }

} // namespace casadi

