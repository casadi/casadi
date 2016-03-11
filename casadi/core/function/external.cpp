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


#include "external_impl.hpp"
#include "../std_vector_tools.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

namespace casadi {
  using namespace std;

  Function external(const string& name, Library li,
                    const Dict& opts) {
    Function ret;
    if (li.has(name + "_simple")) {
      // Simplified, lower overhead external
      ret.assignNode(new SimplifiedExternal(name, li));
    } else {
      // Full information external
      ret.assignNode(new GenericExternal(name, li));
    }
    ret->construct(opts);
    return ret;
  }

  Function external(const string& name, const Dict& opts) {
    return external(name, "./" + name + ".so", opts);
  }

  Function external(const string& name, const string& bin_name,
                    const Dict& opts) {
    return external(name, Library(bin_name), opts);
  }

  Function external(const string& name, const Compiler& compiler,
                    const Dict& opts) {
    return external(name, Library(compiler), opts);
  }

  External::External(const std::string& name, const Library& li)
    : FunctionInternal(name), li_(li) {

    // Increase/decrease reference counter
    incref_ = (signal_t)li_.get(name_ + "_incref");
    decref_ = (signal_t)li_.get(name_ + "_decref");

    // Get number of inputs and outputs
    n_in_ = (getint_t)li_.get(name + "_n_in");
    n_out_ = (getint_t)li_.get(name + "_n_out");

    // Work vector sizes
    work_ = (work_t)li_.get(name_ + "_work");

    // Increase reference counter - external function memory initialized at this point
    if (incref_) incref_();
  }

  SimplifiedExternal::SimplifiedExternal(const std::string& name, const Library& li)
    : External(name, li) {

    // Function for numerical evaluation
    simple_ = (simple_t)li_.get(name_ + "_simple");
  }

  GenericExternal::GenericExternal(const std::string& name, const Library& li)
    : External(name, li) {

    // Functions for retrieving sparsities of inputs and outputs
    sparsity_in_ = (sparsity_t)li_.get(name + "_sparsity_in");
    sparsity_out_ = (sparsity_t)li_.get(name + "_sparsity_out");

    // Function for numerical evaluation
    eval_ = (eval_t)li_.get(name_);
  }

  External::~External() {
    if (decref_) decref_();
    clear_memory();
  }

  size_t External::get_n_in() const {
    if (n_in_) {
      return n_in_();
    } else if (li_.meta().has(name_ + "_N_IN")) {
      return li_.meta().to_int(name_ + "_N_IN");
    } else {
      // Fall back to base class
      return FunctionInternal::get_n_in();
    }
  }

  size_t External::get_n_out() const {
    if (n_out_) {
      return n_out_();
    } else if (li_.meta().has(name_ + "_N_OUT")) {
      return li_.meta().to_int(name_ + "_N_OUT");
    } else {
      // Fall back to base class
      return FunctionInternal::get_n_out();
    }
  }

  Sparsity GenericExternal::get_sparsity_in(int ind) const {
    // Use sparsity retrieval function, if present
    if (sparsity_in_) {
      return Sparsity::compressed(sparsity_in_(ind));
    } else if (li_.meta().has(name_ + "_SPARSITY_IN", ind)) {
      return Sparsity::compressed(li_.meta().to_vector<int>(name_ + "_SPARSITY_IN", ind));
    } else {
      // Fall back to base class
      return FunctionInternal::get_sparsity_in(ind);
    }
  }

  Sparsity GenericExternal::get_sparsity_out(int ind) const {
    // Use sparsity retrieval function, if present
    if (sparsity_out_) {
      return Sparsity::compressed(sparsity_out_(ind));
    } else if (li_.meta().has(name_ + "_SPARSITY_OUT", ind)) {
      return Sparsity::compressed(li_.meta().to_vector<int>(name_ + "_SPARSITY_OUT", ind));
    } else {
      // Fall back to base class
      return FunctionInternal::get_sparsity_out(ind);
    }
  }

  void External::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Reference counting?
    has_refcount_ = li_.has(name_ + "_incref");
    casadi_assert_message(has_refcount_==li_.has(name_ + "_decref"),
                          "External functions must provide functions for both increasing "
                          "and decreasing the reference count, or neither.");

    // Allocate work vectors
    int sz_arg=0, sz_res=0, sz_iw=0, sz_w=0;
    if (work_) {
      int flag = work_(&sz_arg, &sz_res, &sz_iw, &sz_w);
      casadi_assert_message(flag==0, "External: \"work\" failed");
    } else if (li_.meta().has(name_ + "_WORK")) {
      vector<int> v = li_.meta().to_vector<int>(name_ + "_WORK");
      casadi_assert(v.size()==4);
      sz_arg = v[0];
      sz_res = v[1];
      sz_iw = v[2];
      sz_w = v[3];
    }
    alloc_arg(sz_arg);
    alloc_res(sz_res);
    alloc_iw(sz_iw);
    alloc_w(sz_w);
  }

  void SimplifiedExternal::init(const Dict& opts) {
    // Call recursively
    External::init(opts);

    // Arrays for holding inputs and outputs
    alloc_w(n_in() + n_out());
  }

  void GenericExternal::init(const Dict& opts) {
    // Call recursively
    External::init(opts);
  }

  void External::addDependency(CodeGenerator& g) const {
    g.addExternal(signature(name_) + ";");
    if (has_refcount_) {
      g.addExternal("void " + name_ + "_incref(void);");
      g.addExternal("void " + name_ + "_decref(void);");
    }
  }

  std::string External::codegen_name(const CodeGenerator& g) const {
    return name_;
  }

  bool External::hasFullJacobian() const {
    if (FunctionInternal::hasFullJacobian()) return true;
    return li_.has(name_ + "_jac");
  }

  Function External
  ::getFullJacobian(const std::string& name, const Dict& opts) {
    if (hasFullJacobian()) {
      return external(name, li_, opts);
    } else {
      return FunctionInternal::getFullJacobian(name, opts);
    }
  }

  Function External
  ::get_forward(const std::string& name, int nfwd, Dict& opts) {
    // Consistency check
    int n=1;
    while (n<nfwd) n*=2;
    if (n!=nfwd || nfwd>get_n_forward()) {
      casadi_error("Internal error: Refactoring needed, cf. #1055");
    }
    return external(name, li_, opts);
  }

  int External::get_n_forward() const {
    // Will try 64, 32, 16, 8, 4, 2, 1 directions
    for (int i=64; i>0; i/=2) {
      stringstream ss;
      ss << "fwd" << i << "_" << name_;
      if (li_.has(ss.str())) return i;
    }
    return 0;
  }

  Function External
  ::get_reverse(const std::string& name, int nadj, Dict& opts) {
    // Consistency check
    int n=1;
    while (n<nadj) n*=2;
    if (n!=nadj || nadj>get_n_reverse()) {
      casadi_error("Internal error: Refactoring needed, cf. #1055");
    }
    return external(name, li_, opts);
  }

  int External::get_n_reverse() const {
    // Will try 64, 32, 16, 8, 4, 2, 1 directions
    for (int i=64; i>0; i/=2) {
      stringstream ss;
      ss << "adj" << i << "_" << name_;
      if (li_.has(ss.str())) return i;
    }
    return 0;
  }

} // namespace casadi

