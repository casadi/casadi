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


#include "compiler.hpp"
#include "compiler_internal.hpp"

using namespace std;
namespace casadi {

  Compiler::Compiler() {
  }

  Compiler::Compiler(const std::string& name,
                           const std::string& compiler,
                           const Dict& opts) {
    assignNode(CompilerInternal::getPlugin(compiler).creator(name));
    setOption(opts);
    init();
  }

  CompilerInternal* Compiler::operator->() {
    return static_cast<CompilerInternal*>(OptionsFunctionality::operator->());
  }

  const CompilerInternal* Compiler::operator->() const {
    return static_cast<const CompilerInternal*>(OptionsFunctionality::operator->());
  }

  bool Compiler::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const CompilerInternal*>(ptr)!=0;
  }

  bool Compiler::hasPlugin(const std::string& name) {
    return CompilerInternal::hasPlugin(name);
  }

  void Compiler::loadPlugin(const std::string& name) {
    CompilerInternal::loadPlugin(name);
  }

  std::string Compiler::doc(const std::string& name) {
    return CompilerInternal::getPlugin(name).doc;
  }

  std::string Compiler::plugin_name() const {
    return (*this)->plugin_name();
  }

  void* Compiler::getFunction(const std::string& symname) {
    return (*this)->getFunction(symname);
  }

} // namespace casadi
