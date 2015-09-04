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


#include "jit_function.hpp"
#include "jit_function_internal.hpp"
#include "sx_function.hpp"
#include "mx_function.hpp"

using namespace std;
namespace casadi {

  JitFunction::JitFunction() {
  }

  JitFunction::JitFunction(
      const std::string& compiler,
      const Function& f,
      const Dict& opts) {
    assignNode(JitFunctionInternal::instantiatePlugin(compiler, f));
    setOption("name", f.getOption("name"));
    setOption(opts);
    init();
  }

  JitFunctionInternal* JitFunction::operator->() {
    return static_cast<JitFunctionInternal*>(Function::operator->());
  }

  const JitFunctionInternal* JitFunction::operator->() const {
    return static_cast<const JitFunctionInternal*>(Function::operator->());
  }

  bool JitFunction::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const JitFunctionInternal*>(ptr)!=0;
  }

  bool JitFunction::hasPlugin(const std::string& name) {
    return JitFunctionInternal::hasPlugin(name);
  }

  void JitFunction::loadPlugin(const std::string& name) {
    JitFunctionInternal::loadPlugin(name);
  }

  std::string JitFunction::doc(const std::string& name) {
    return JitFunctionInternal::getPlugin(name).doc;
  }

} // namespace casadi
