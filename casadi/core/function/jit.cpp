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


#include "jit.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

namespace casadi {
  using namespace std;

  Jit::Jit(const std::string& name, int n_in, int n_out,
           const std::string& body, const Dict& opts)
  : FunctionInternal(name), n_in_(n_in), n_out_(n_out), body_(body) {
    addOption("jac", OT_STRING, GenericType(), "Function body for Jacobian");
    addOption("hess", OT_STRING, GenericType(), "Function body for Hessian");

    // Jit by default
    setOption("jit", true);

    // Arrays for holding inputs and outputs
    alloc_w(n_in + n_out);
  }

  void Jit::init() {
    // Call the initialization method of the base class
    FunctionInternal::init();

    // Read options
    if (hasSetOption("jac")) jac_body_ = option("jac").toString();
    if (hasSetOption("hess")) hess_body_ = option("hess").toString();
  }

  Jit::~Jit() {
  }

  void Jit::generateBody(CodeGenerator& g) const {
    g.body << body_;
  }

  bool Jit::hasFullJacobian() const {
    return !jac_body_.empty();
  }

  Function Jit::getFullJacobian(const std::string& name, const Dict& opts) {
    // Create a JIT-function for the Jacobian
    Dict jit_opts;
    if (!hess_body_.empty()) jit_opts["jac"] = hess_body_;
    Function fcn = jit(name + "_jit", n_in_, n_in_*n_out_, jac_body_, jit_opts);

    // Wrap in an MX function
    std::vector<MX> arg = fcn.mx_in();
    std::vector<MX> res = fcn(arg);
    MX J = reshape(vertcat(res), n_in_, n_out_);
    return Function(name, arg, {J}, opts);
  }

} // namespace casadi
