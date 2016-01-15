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


#include "compiler_internal.hpp"

using namespace std;
namespace casadi {

  CompilerInternal::CompilerInternal(const std::string& name) : name_(name) {
  }

  CompilerInternal::~CompilerInternal() {
  }

  void CompilerInternal::print(ostream &stream) const {
    stream << "Compiler" << endl;
  }

  void CompilerInternal::repr(ostream &stream) const {
    stream << "Compiler" << endl;
  }

  std::map<std::string, CompilerInternal::Plugin> CompilerInternal::solvers_;

  const std::string CompilerInternal::infix_ = "compiler";

  Options CompilerInternal::options_
  = {{},
     {{}
     }
  };

  void CompilerInternal::construct(const Dict& opts) {
    // Get a reference to the options structure
    const Options& options = get_options();

    // Make sure all options exist and have the correct type
    for (auto&& op : opts) {
      const Options::Entry* entry = options.find(op.first);
      casadi_assert_message(entry!=0, "No such option: " + op.first);
      casadi_assert_message(op.second.can_cast_to(entry->type),
                            "Illegal type for " + op.first);
    }

    init(opts);
  }

} // namespace casadi
