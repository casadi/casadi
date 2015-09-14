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


#include "jit_compiler.hpp"
#include "jit_compiler_internal.hpp"
#include "sx_function.hpp"
#include "mx_function.hpp"

using namespace std;
namespace casadi {

  JitCompiler::JitCompiler() {
  }

  JitCompiler::JitCompiler(
      const std::string& compiler,
      const Function& f,
      const Dict& opts) {
    assignNode(JitCompilerInternal::instantiatePlugin(compiler, f));
    setOption("name", f.getOption("name"));
    setOption(opts);
    init();
  }

  JitCompilerInternal* JitCompiler::operator->() {
    return static_cast<JitCompilerInternal*>(Function::operator->());
  }

  const JitCompilerInternal* JitCompiler::operator->() const {
    return static_cast<const JitCompilerInternal*>(Function::operator->());
  }

  bool JitCompiler::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const JitCompilerInternal*>(ptr)!=0;
  }

  bool JitCompiler::hasPlugin(const std::string& name) {
    return JitCompilerInternal::hasPlugin(name);
  }

  void JitCompiler::loadPlugin(const std::string& name) {
    JitCompilerInternal::loadPlugin(name);
  }

  std::string JitCompiler::doc(const std::string& name) {
    return JitCompilerInternal::getPlugin(name).doc;
  }

  std::string JitCompiler::plugin_name() const {
    return (*this)->plugin_name();
  }

  void* JitCompiler::getFunction(const std::string& symname) {
    return (*this)->getFunction(symname);
  }

} // namespace casadi
