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


#include "jit_impl.hpp"
#include "casadi_file.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

namespace casadi {
  using namespace std;

  Function jit(const std::string& name, int n_in, int n_out,
               const std::string& body, const Dict& opts) {
    return Function::create(new Jit(name, n_in, n_out, body, opts), opts);
  }

  Function jit(const ParsedFile& file) {
    // Function name
    string name = file.to_string("NAME");

    // Number of inputs (default 1)
    int n_in = file.has("N_IN") ? file.to_int("N_IN") : 1;

    // Number of outputs (default 1)
    int n_out = file.has("N_OUT") ? file.to_int("N_OUT") : 1;

    // Function body
    string body = file.to_text("BODY");

    // Read options
    Dict opts;
    if (file.has("JAC_BODY")) opts["jac"] = file.to_text("JAC_BODY");
    if (file.has("HESS_BODY")) opts["hess"] = file.to_text("HESS_BODY");
    if (file.has("JIT")) opts["jit"] = file.to_int("JIT");

    // Create function object
    return jit(name, n_in, n_out, body, opts);
  }

  Jit::Jit(const std::string& name, int n_in, int n_out,
           const std::string& body, const Dict& opts)
    : FunctionInternal(name), n_in_(n_in), n_out_(n_out), body_(body) {

    // Default options
    jit_ = true; // override default
  }

  Options Jit::options_
  = {{&FunctionInternal::options_},
     {{"jac",
       {OT_STRING,
        "Function body for Jacobian"}},
      {"hess",
       {OT_STRING,
        "Function body for Hessian"}}
     }
  };

  void Jit::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Read options
    for (auto&& op : opts) {
      if (op.first=="jac") {
        jac_body_ = op.second.to_string();
      } else if (op.first=="hess") {
        hess_body_ = op.second.to_string();
      }
    }

    // Arrays for holding inputs and outputs
    alloc_w(n_in() + n_out());
  }

  Jit::~Jit() {
  }

  void Jit::codegen_body(CodeGenerator& g) const {
    g << body_;
  }

  bool Jit::has_jacobian() const {
    return !jac_body_.empty();
  }

  Function Jit::get_jacobian(const std::string& name,
                                const std::vector<std::string>& inames,
                                const std::vector<std::string>& onames,
                                const Dict& opts) const {
    // Create a JIT-function for the Jacobian
    Dict jit_opts;
    if (!hess_body_.empty()) jit_opts["jac"] = hess_body_;
    Function fcn = jit(name + "_jit", n_in_ + n_out_, n_in_*n_out_, jac_body_, jit_opts);

    // Wrap in an MX function FIXME(@jaeandersson)
    std::vector<MX> arg = fcn.mx_in();
    std::vector<MX> res = fcn(arg);
    MX J = reshape(vertcat(res), n_in_, n_out_);
    return Function(name, arg, {J}, opts);
  }

} // namespace casadi
