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

  TranslatorInternal::TranslatorInternal(const std::string& name)
    : name_(name), verbose_(false) {
  }

  TranslatorInternal::~TranslatorInternal() {
  }

  const Options TranslatorInternal::options_
  = {{},
     {{"verbose",
       {OT_BOOL,
        "Verbose evaluation -- for debugging"}}}
  };

  void TranslatorInternal::construct(const Dict& opts) {
    // Read options
    for (auto&& op : opts) {
      if (op.first=="verbose") {
        verbose_ = op.second;
      }
    }
    init(opts);
  }

  void TranslatorInternal::init(const Dict& opts) {
    // Default implementation does nothing
  }

  void TranslatorInternal::disp(std::ostream& stream, bool more) const {
    stream << "Translator (";
    stream << "plugin: " << plugin_name();
    stream << ", name: " << name_;
    stream << ")";
  }

  std::map<std::string, TranslatorInternal::Plugin> TranslatorInternal::solvers_;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
  std::mutex TranslatorInternal::mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

  const std::string TranslatorInternal::infix_ = "translator";

} // namespace casadi
