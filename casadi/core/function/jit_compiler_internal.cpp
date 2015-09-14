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


#include "jit_compiler_internal.hpp"
#include "mx_function.hpp"
#include "sx_function.hpp"

using namespace std;
namespace casadi {

  JitCompilerInternal::JitCompilerInternal(const Function& f) : f_(f) {

    // Set default options
    setOption("name", "unnamed JIT function"); // name of the function
    if (f.hasSetOption("input_scheme")) {
      setOption("input_scheme", f.getOption("input_scheme"));
    }
    if (f.hasSetOption("output_scheme")) {
      setOption("output_scheme", f.getOption("output_scheme"));
    }
  }

  JitCompilerInternal::~JitCompilerInternal() {

  }

  std::map<std::string, JitCompilerInternal::Plugin> JitCompilerInternal::solvers_;

  const std::string JitCompilerInternal::infix_ = "jitcompiler";

  void JitCompilerInternal::init() {
    // Initialize the base classes
    FunctionInternal::init();
  }

} // namespace casadi
