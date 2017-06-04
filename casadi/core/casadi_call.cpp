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


#include "casadi_call.hpp"
#include "function_internal.hpp"
#include "std_vector_tools.hpp"

using namespace std;

namespace casadi {

  MX Call::projectArg(const MX& x, const Sparsity& sp, int i) {
    if (x.size()==sp.size()) {
      // Insert sparsity projection nodes if needed
      return project(x, sp);
    } else {
      // Different dimensions
      if (x.is_empty() || sp.is_empty()) { // NOTE: To permissive?
        // Replace nulls with zeros of the right dimension
        return MX::zeros(sp);
      } else if (x.is_scalar()) {
        // Scalar argument means set all
        return MX(sp, x);
      } else if (x.size1()==sp.size2() && x.size2()==sp.size1() && sp.is_vector()) {
        // Transposed vector
        return projectArg(x.T(), sp, i);
      } else {
        // Mismatching dimensions
        casadi_error("Cannot create function call node: Dimension mismatch for argument "
                     << i << ". Argument has shape " << x.size()
                     << " but function input has shape " << sp.size());
      }
    }
  }

  Call::Call(const Function& fcn, const vector<MX>& arg) : fcn_(fcn) {

    // Number inputs and outputs
    int num_in = fcn.n_in();
    casadi_assert_message(arg.size()==num_in, "Argument list length (" << arg.size()
                          << ") does not match number of inputs (" << num_in << ")"
                          << " for function " << fcn.name());

    // Create arguments of the right dimensions and sparsity
    vector<MX> arg1(num_in);
    for (int i=0; i<num_in; ++i) {
      arg1[i] = projectArg(arg[i], fcn_.sparsity_in(i), i);
    }
    set_dep(arg1);
    set_sparsity(Sparsity::scalar());
  }

  std::string Call::print(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << fcn_.name() << "(";
    for (int i=0; i<ndep(); ++i) {
      if (i!=0) ss << ", ";
      ss << arg.at(i);
    }
    ss << ")";
    return ss.str();
  }

  void Call::eval(const double** arg, double** res, int* iw, double* w, int mem) const {
    fcn_(arg, res, iw, w, mem);
  }

  int Call::nout() const {
    return fcn_.n_out();
  }

  const Sparsity& Call::sparsity(int oind) const {
    return fcn_.sparsity_out(oind);
  }

  void Call::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) const {
    fcn_(arg, res, iw, w, mem);
  }

  void Call::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res = create(fcn_, arg);
  }

  void Call::ad_forward(const vector<vector<MX> >& fseed,
                     vector<vector<MX> >& fsens) const {
    // Nondifferentiated inputs and outputs
    vector<MX> arg(ndep());
    for (int i=0; i<arg.size(); ++i) arg[i] = dep(i);
    vector<MX> res(nout());
    for (int i=0; i<res.size(); ++i) res[i] = get_output(i);

    // Call the cached functions
    fcn_->call_forward(arg, res, fseed, fsens, false, false);
  }

  void Call::ad_reverse(const vector<vector<MX> >& aseed,
                     vector<vector<MX> >& asens) const {
    // Nondifferentiated inputs and outputs
    vector<MX> arg(ndep());
    for (int i=0; i<arg.size(); ++i) arg[i] = dep(i);
    vector<MX> res(nout());
    for (int i=0; i<res.size(); ++i) res[i] = get_output(i);

    // Call the cached functions
    vector<vector<MX> > v;
    fcn_->call_reverse(arg, res, aseed, v, false, false);
    for (int i=0; i<v.size(); ++i) {
      for (int j=0; j<v[i].size(); ++j) {
        if (!v[i][j].is_empty()) { // TODO(@jaeandersson): Hack
          asens[i][j] += v[i][j];
        }
      }
    }
  }

  void Call::sp_forward(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const {
    fcn_(arg, res, iw, w, mem);
  }

  void Call::sp_reverse(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const {
    fcn_.rev(arg, res, iw, w, mem);
  }

  void Call::addDependency(CodeGenerator& g) const {
    fcn_->addDependency(g);
  }

  bool Call::has_refcount() const {
    return fcn_->has_refcount_;
  }

  void Call::generate(CodeGenerator& g, const std::string& mem,
                      const vector<int>& arg, const vector<int>& res) const {
    if (fcn_->simplifiedCall()) {

      // Collect input arguments
      for (int i=0; i<arg.size(); ++i) {
        g << "w[" << i << "]=" << g.workel(arg[i]) << ";\n";
      }

      // Call function
      g << g(fcn_, "w", "w+"+g.to_string(arg.size())) << ";\n";

      // Collect output arguments
      for (int i=0; i<res.size(); ++i) {
        if (res[i]>=0) {
          g << g.workel(res[i]) << "=w[" << (arg.size()+i) << "];\n";
        }
      }
    } else {
      // Collect input arguments
      g.local("arg1", "const real_t", "**");
      for (int i=0; i<arg.size(); ++i) {
        g << "arg1[" << i << "]=" << g.work(arg[i], fcn_.nnz_in(i)) << ";\n";
      }

      // Collect output arguments
      g.local("res1", "real_t", "**");
      for (int i=0; i<res.size(); ++i) {
        g << "res1[" << i << "]=" << g.work(res[i], fcn_.nnz_out(i)) << ";\n";
      }

      // Call function
      g << "if (" << g(fcn_, "arg1", "res1", "iw", "w") << ") return 1;\n";
    }
  }

  void Call::codegen_incref(CodeGenerator& g, std::set<void*>& added) const {
    if (has_refcount()) {
      auto i = added.insert(fcn_.get());
      if (i.second) { // prevent duplicate calls
        g << fcn_->codegen_name(g) << "_incref();\n";
      }
    }
  }

  void Call::codegen_decref(CodeGenerator& g, std::set<void*>& added) const {
    if (has_refcount()) {
      auto i = added.insert(fcn_.get());
      if (i.second) { // prevent duplicate calls
        g << fcn_->codegen_name(g) << "_decref();\n";
      }
    }
  }

  size_t Call::sz_arg() const {
    return fcn_.sz_arg();
  }

  size_t Call::sz_res() const {
    return fcn_.sz_res();
  }

  size_t Call::sz_iw() const {
    return fcn_.sz_iw();
  }

  size_t Call::sz_w() const {
    return fcn_.sz_w();
  }

  std::vector<MX> Call::create(const Function& fcn, const std::vector<MX>& arg) {
    return MX::createMultipleOutput(new Call(fcn, arg));
  }

} // namespace casadi
