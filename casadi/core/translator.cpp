/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#include "translator_internal.hpp"

namespace casadi {

  Translator::Translator() {
  }

  Translator::Translator(const std::string& name, const Dict& opts) {
    own(TranslatorInternal::getPlugin(name).creator());
    (*this)->construct(opts);
  }

  TranslatorInternal* Translator::operator->() {
    return static_cast<TranslatorInternal*>(SharedObject::operator->());
  }

  const TranslatorInternal* Translator::operator->() const {
    return static_cast<const TranslatorInternal*>(SharedObject::operator->());
  }

  bool Translator::test_cast(const SharedObjectInternal* ptr) {
    return dynamic_cast<const TranslatorInternal*>(ptr)!=nullptr;
  }

  bool Translator::has_plugin(const std::string& name) {
    return TranslatorInternal::has_plugin(name);
  }

  void Translator::load_plugin(const std::string& name) {
    TranslatorInternal::load_plugin(name);
  }

  std::string Translator::doc(const std::string& name) {
    return TranslatorInternal::getPlugin(name).doc;
  }

  std::string Translator::plugin_name() const {
    return (*this)->plugin_name();
  }

  void Translator::load(const std::string& filename) {
    (*this)->load(filename);
  }

  void Translator::load(const Function& f) {
    (*this)->load(f);
  }

  void Translator::set_dimension(const std::string& name, casadi_int dim) {
    (*this)->set_dimension(name, dim);
  }

  Function Translator::create(const std::string& name) {
    return (*this)->create(name);
  }

  void Translator::save(const std::string& filename) {
    (*this)->save(filename);
  }

  Translator Translator::create(TranslatorInternal* node) {
    Translator ret;
    ret.own(node);
    return ret;
  }

  Translator Translator::create(TranslatorInternal* node, const Dict& opts) {
    Translator ret = Translator::create(node);
    ret->construct(opts);
    return ret;
  }

  bool has_translator(const std::string& name) {
    return Translator::has_plugin(name);
  }

  void load_translator(const std::string& name) {
    Translator::load_plugin(name);
  }

  std::string doc_translator(const std::string& name) {
    return Translator::doc(name);
  }

} // namespace casadi
