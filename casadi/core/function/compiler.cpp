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
    if (compiler=="none") {
      assignNode(new CompilerInternal(name));
    } else if (compiler=="dll") {
      assignNode(new DllLibrary(name));
    } else {
      assignNode(CompilerInternal::getPlugin(compiler).creator(name));
    }
    (*this)->construct(opts);
  }

  CompilerInternal* Compiler::operator->() {
    return static_cast<CompilerInternal*>(SharedObject::operator->());
  }

  const CompilerInternal* Compiler::operator->() const {
    return static_cast<const CompilerInternal*>(SharedObject::operator->());
  }

  bool Compiler::test_cast(const SharedObjectNode* ptr) {
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

  bool Compiler::has_function(const std::string& symname) const {
    return (*this)->has_function(symname);
  }

  signal_t Compiler::get_function(const std::string& symname) {
    return (*this)->get_function(symname);
  }

  bool Library::has_meta(const std::string& cmd, int ind) const {
    return (*this)->meta().has(cmd, ind);
  }

  std::string Library::get_meta(const std::string& cmd, int ind) const {
    return (*this)->meta().to_text(cmd, ind);
  }

  Library::Library() {
  }

  Library::Library(const std::string& bin_name) {
    Compiler compiler(bin_name, "dll");
    assignNode(new JitLibrary(compiler));
  }

  Library::Library(const Compiler& compiler) {
    assignNode(new JitLibrary(compiler));
  }

  LibraryInternal* Library::operator->() {
    return static_cast<LibraryInternal*>(SharedObject::operator->());
  }

  const LibraryInternal* Library::operator->() const {
    return static_cast<const LibraryInternal*>(SharedObject::operator->());
  }

  bool Library::has(const std::string& sym) const {
    return (*this)->has(sym);
  }

  signal_t Library::get(const std::string& sym) {
    return (*this)->get(sym);
  }

} // namespace casadi
