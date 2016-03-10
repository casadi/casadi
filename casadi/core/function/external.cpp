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
    return external("./" + name + ".so", name, opts);
  }

  Function external(const string& name, const string& bin_name,
                    const Dict& opts) {
    return external(name, Library(bin_name), opts);
  }

  Function external(const string& name, const Compiler& compiler,
                    const Dict& opts) {
    return external(name, Library(compiler), opts);
  }

  Options External::options_
  = {{&FunctionInternal::options_},
     {{"int_data",
       {OT_INTVECTOR,
        "Integer data vector to be passed to the external function"}},
      {"real_data",
       {OT_DOUBLEVECTOR,
        "Real data vector to be passed to the external function"}}
     }
  };

  void External::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    /** \brief Parameter lengths for consistency check */
    int n_int = int_data_.size();
    int n_real = real_data_.size();

    // Read options
    for (auto&& op : opts) {
      if (op.first=="int_data") {
        int_data_ = op.second;
      } else if (op.first=="real_data") {
        real_data_ = op.second;
      }
    }

    // Consistency check
    casadi_assert_message(int_data_.size()==n_int, "'int_data' has wrong length");
    casadi_assert_message(real_data_.size()==n_real, "'real_data' has wrong length");

    // Get number of temporaries
    int sz_arg, sz_res, sz_iw, sz_w;
    work_t work = (work_t)li_.get(name_ + "_work");
    if (work!=0) {
      int flag = work(&sz_arg, &sz_res, &sz_iw, &sz_w);
      casadi_assert_message(flag==0, "External: \"work\" failed");
    } else {
      // No work vectors
      sz_arg = n_in();
      sz_res = n_out();
      sz_iw = sz_w = 0;
    }

    // Function for allocating memory, if any
    allocmem_ = (allocmem_t)li_.get(name_ + "_alloc");

    // Function for freeing memory, if any
    freemem_ = (freemem_t)li_.get(name_ + "_free");

    alloc_arg(sz_arg);
    alloc_res(sz_res);
    alloc_iw(sz_iw);
    alloc_w(sz_w);
  }

  SimplifiedExternal::SimplifiedExternal(const std::string& name, const Library& li)
    : External(name, li) {
  }

  void SimplifiedExternal::init(const Dict& opts) {
    // Call recursively
    External::init(opts);

    // Function for numerical evaluation
    string eval_s = this->name_ + "_simple";
    eval_ = (simple_t)li_.get(eval_s);
    casadi_assert_message(eval_!=0, "SimplifiedExternal: no \"" + eval_s + "\" found");

    // Arrays for holding inputs and outputs
    this->alloc_w(this->n_in() + this->n_out());
  }

  size_t External::get_n_in() const {
    if (init_) {
      int n_in;
      int flag = init_(&n_in, 0, 0, 0);
      casadi_assert_message(flag==0, "External: \"init\" failed");
      return n_in;
    } else if (li_.meta().has(name_ + "_N_IN")) {
      return li_.meta().to_int(name_ + "_N_IN");
    } else {
      // Fall back to base class
      return External::get_n_in();
    }
  }

  size_t External::get_n_out() const {
    if (init_) {
      int n_out;
      int flag = init_(0, &n_out, 0, 0);
      casadi_assert_message(flag==0, "External: \"init\" failed");
      return n_out;
    } else if (li_.meta().has(name_ + "_N_OUT")) {
      return li_.meta().to_int(name_ + "_N_OUT");
    } else {
      // Fall back to base class
      return External::get_n_out();
    }
  }

  Sparsity External::get_sparsity_in(int ind) const {
    // Use sparsity retrieval function, if present
    if (sparsity_) {
      // Get sparsity pattern in CCS format
      int nrow, ncol;
      const int *colind, *row;
      int flag = sparsity_(ind, &nrow, &ncol, &colind, &row);
      casadi_assert_message(flag==0, "External: \"sparsity\" failed");
      return Sparsity(nrow, ncol, colind, row);
    } else if (li_.meta().has(name_ + "_SPARSITY_IN", ind)) {
      const ParsedFile& m = li_.meta();
      vector<int> v = m.to_vector<int>(name_ + "_SPARSITY_IN", ind);
      return Sparsity::compressed(v);
    } else {
      // Fall back to base class
      return External::get_sparsity_in(ind);
    }
  }

  Sparsity External::get_sparsity_out(int ind) const {
    // Use sparsity retrieval function, if present
    if (sparsity_) {
      // Get sparsity pattern in CCS format
      int nrow, ncol;
      const int *colind, *row;
      int flag = sparsity_(ind+get_n_in(), &nrow, &ncol, &colind, &row);
      casadi_assert_message(flag==0, "External: \"sparsity\" failed");
      return Sparsity(nrow, ncol, colind, row);
    } else if (li_.meta().has(name_ + "_SPARSITY_OUT", ind)) {
      const ParsedFile& m = li_.meta();
      vector<int> v = m.to_vector<int>(name_ + "_SPARSITY_OUT", ind);
      return Sparsity::compressed(v);
    } else {
      // Fall back to base class
      return External::get_sparsity_out(ind);
    }
  }

  GenericExternal::GenericExternal(const std::string& name, const Library& li)
    : External(name, li) {

    // Function for retrieving sparsities of inputs and outputs
    sparsity_ = (sparsity_t)li_.get(name + "_sparsity");
  }

  External::External(const std::string& name, const Library& li)
    : FunctionInternal(name), li_(li) {

    // Initialize and get the number of inputs and outputs
    // Function for creating an object instance
    string init_s = name + "_init";
    int n_int=0, n_real=0;
    init_ = (init_t)li_.get(init_s);
    if (init_!=0) {
      int flag = init_(0, 0, &n_int, &n_real);
      casadi_assert_message(flag==0, "External: \"init\" failed");
    } else {
      // Fall back to reading meta information
      if (li_.meta().has(name + "_N_INT")) n_int=li_.meta().to_int(name + "_N_INT");
      if (li_.meta().has(name + "_N_REAL")) n_real=li_.meta().to_int(name + "_N_REAL");
    }
    int_data_.resize(n_int);
    real_data_.resize(n_real);
  }

  External::~External() {
    if (freemem_) {
      freemem_(0);
    }
  }

  void GenericExternal::init(const Dict& opts) {
    // Call recursively
    External::init(opts);

    // Function for numerical evaluation
    eval_ = (eval_t)li_.get(name_);
  }

  void SimplifiedExternal::
  simple(const double* arg, double* res) {
    casadi_assert_message(eval_!=0, "Numerical evaluation not possible");
    eval_(arg, res);
  }

  void GenericExternal::
  eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
    casadi_assert_message(eval_!=0, "Numerical evaluation not possible");
    int flag = eval_(mem, arg, res, iw, w);
    if (flag) throw CasadiException("External: \""+this->name_+"\" failed");
  }

  void External::addDependency(CodeGenerator& g) const {
    g.addExternal(signature(this->name_) + ";");
  }

  std::string SimplifiedExternal::
  simple_call(const CodeGenerator& g,
              const std::string& arg, const std::string& res) const {
    // Create a function call
    stringstream ss;
    ss << this->name_ << "(" << arg << ", " << res << ")";
    return ss.str();
  }

  std::string GenericExternal::
  generic_call(const CodeGenerator& g, const std::string& mem,
               const std::string& arg, const std::string& res,
               const std::string& iw, const std::string& w) const {

    // Create a function call
    stringstream ss;
    ss << this->name_ << "(" << mem << ", " << arg << ", " << res << ", "
       << iw << ", " << w << ")";
    return ss.str();
  }

  bool External::hasFullJacobian() const {
    if (FunctionInternal::hasFullJacobian()) return true;
    return const_cast<Library&>(li_).has(name_ + "_jac");
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
    init_t fwd_init;
    for (int i=64; i>0; i/=2) {
      stringstream ss;
      ss << "fwd" << i << "_" << name_;
      fwd_init = (init_t)const_cast<Library&>(li_).get(ss.str());
      if (fwd_init!=0) return i;
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
    init_t adj_init;
    for (int i=64; i>0; i/=2) {
      stringstream ss;
      ss << "adj" << i << "_" << name_;
      adj_init = (init_t)const_cast<Library&>(li_).get(ss.str());
      if (adj_init!=0) return i;
    }
    return 0;
  }

} // namespace casadi

