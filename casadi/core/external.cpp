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


#include "external_impl.hpp"
#include "casadi_misc.hpp"
#include "serializing_stream.hpp"
#include <casadi/config.h>

#include <fstream>
#include <iostream>
#include <sstream>

// Set default shared library suffix
#ifndef SHARED_LIBRARY_SUFFIX
#define SHARED_LIBRARY_SUFFIX CASADI_SHARED_LIBRARY_SUFFIX
#endif // SHARED_LIBRARY_SUFFIX

namespace casadi {

Function external(const std::string& name, const Importer& li,
                  const Dict& opts) {
  return Function::create(new GenericExternal(name, li), opts);
}

Function external(const std::string& name, const Dict& opts) {
  return external(name, "./" + name + SHARED_LIBRARY_SUFFIX, opts);
}

Function external(const std::string& name, const std::string& bin_name,
                  const Dict& opts) {
  return external(name, Importer(bin_name, "dll"), opts);
}

External::External(const std::string& name, const Importer& li)
  : FunctionInternal(name), li_(li) {

  External::init_external();
}

bool External::any_symbol_found() const {
  return incref_ || decref_ || get_default_in_ ||
    get_n_in_ || get_n_out_ || get_name_in_ ||
    get_name_out_ || work_;
}

void External::init_external() {
  // Increasing/decreasing reference counter
  incref_ = (signal_t)li_.get_function(name_ + "_incref");
  decref_ = (signal_t)li_.get_function(name_ + "_decref");

  casadi_assert(static_cast<bool>(incref_) == static_cast<bool>(decref_),
    "External must either define both incref and decref or neither.");

  // Getting default arguments
  get_default_in_ = (default_t)li_.get_function(name_ + "_default_in");

  // Getting number of inputs and outputs
  get_n_in_ = (getint_t)li_.get_function(name_ + "_n_in");
  get_n_out_ = (getint_t)li_.get_function(name_ + "_n_out");

  // Getting names of inputs and outputs
  get_name_in_ = (name_t)li_.get_function(name_ + "_name_in");
  get_name_out_ = (name_t)li_.get_function(name_ + "_name_out");

  // Work vector sizes
  work_ = (work_t)li_.get_function(name_ + "_work");

  // Increase reference counter - external function memory initialized at this point
  if (incref_) incref_();
}

GenericExternal::GenericExternal(const std::string& name, const Importer& li)
  : External(name, li) {

  GenericExternal::init_external();
}

bool GenericExternal::any_symbol_found() const {
  return External::any_symbol_found() ||
    get_sparsity_in_ || get_sparsity_out_ ||
    get_diff_in_ || get_diff_out_ ||
    checkout_ || release_ ||
    eval_;
}

void GenericExternal::init_external() {

  // Functions for retrieving sparsities of inputs and outputs
  get_sparsity_in_ = (sparsity_t)li_.get_function(name_ + "_sparsity_in");
  get_sparsity_out_ = (sparsity_t)li_.get_function(name_ + "_sparsity_out");

  // Differentiability of inputs and outputs
  get_diff_in_ = (diff_t)li_.get_function(name_ + "_diff_in");
  get_diff_out_ = (diff_t)li_.get_function(name_ + "_diff_out");

  // Memory management functions
  checkout_ = (casadi_checkout_t) li_.get_function(name_ + "_checkout");
  release_ = (casadi_release_t) li_.get_function(name_ + "_release");

  casadi_assert(static_cast<bool>(checkout_) == static_cast<bool>(release_),
    "External must either define both checkout and release or neither.");

  // Function for numerical evaluation
  eval_ = (eval_t)li_.get_function(name_);

  // Sparsity patterns of Jacobians
  get_jac_sparsity_ = (sparsity_t)li_.get_function("jac_" + name_ + "_sparsity_out");
}

External::~External() {
  if (decref_) decref_();
  clear_mem();
}

size_t External::get_n_in() {
  if (get_n_in_) {
    return get_n_in_();
  } else if (li_.has_meta(name_ + "_N_IN")) {
    return li_.meta_int(name_ + "_N_IN");
  } else {
    // Fall back to base class
    return FunctionInternal::get_n_in();
  }
}

size_t External::get_n_out() {
  if (get_n_out_) {
    return get_n_out_();
  } else if (li_.has_meta(name_ + "_N_OUT")) {
    return li_.meta_int(name_ + "_N_OUT");
  } else {
    // Fall back to base class
    return FunctionInternal::get_n_out();
  }
}

double External::get_default_in(casadi_int i) const {
  if (get_default_in_) {
    return get_default_in_(i);
  } else {
    // Fall back to base class
    return FunctionInternal::get_default_in(i);
  }
}

std::string External::get_name_in(casadi_int i) {
  if (get_name_in_) {
    // Use function pointer
    const char* n = get_name_in_(i);
    casadi_assert(n!=nullptr, "Error querying input name");
    return n;
  } else if (li_.has_meta(name_ + "_NAME_IN", i)) {
    // Read meta
    return li_.meta_string(name_ + "_NAME_IN", i);
  } else {
    // Default name
    return FunctionInternal::get_name_in(i);
  }
}

std::string External::get_name_out(casadi_int i) {
  if (get_name_out_) {
    // Use function pointer
    const char* n = get_name_out_(i);
    casadi_assert(n!=nullptr, "Error querying output name");
    return n;
  } else if (li_.has_meta(name_ + "_NAME_OUT", i)) {
    // Read meta
    return li_.meta_string(name_ + "_NAME_OUT", i);
  } else {
    // Default name
    return FunctionInternal::get_name_out(i);
  }
}

Sparsity GenericExternal::get_sparsity_in(casadi_int i) {
  // Use sparsity retrieval function, if present
  if (get_sparsity_in_) {
    return Sparsity::compressed(get_sparsity_in_(i));
  } else if (li_.has_meta(name_ + "_SPARSITY_IN", i)) {
    return Sparsity::compressed(li_.meta_vector<casadi_int>(name_ + "_SPARSITY_IN", i));
  } else {
    // Fall back to base class
    return FunctionInternal::get_sparsity_in(i);
  }
}

Sparsity GenericExternal::get_sparsity_out(casadi_int i) {
  // Use sparsity retrieval function, if present
  if (get_sparsity_out_) {
    return Sparsity::compressed(get_sparsity_out_(i));
  } else if (li_.has_meta(name_ + "_SPARSITY_OUT", i)) {
    return Sparsity::compressed(li_.meta_vector<casadi_int>(name_ + "_SPARSITY_OUT", i));
  } else {
    // Fall back to base class
    return FunctionInternal::get_sparsity_out(i);
  }
}

bool GenericExternal::has_jac_sparsity(casadi_int oind, casadi_int iind) const {
  // Flat index
  casadi_int ind = iind + oind * n_in_;
  // Jacobian sparsity pattern known?
  if (get_jac_sparsity_ || li_.has_meta("JAC_" + name_ + "_SPARSITY_OUT", ind)) {
    return true;
  } else {
    // Fall back to base class
    return FunctionInternal::has_jac_sparsity(oind, iind);
  }
}

Sparsity GenericExternal::get_jac_sparsity(casadi_int oind, casadi_int iind,
    bool symmetric) const {
  // Flat index
  casadi_int ind = iind + oind * n_in_;
  // Use sparsity retrieval function, if present
  if (get_jac_sparsity_) {
    return Sparsity::compressed(get_jac_sparsity_(ind));
  } else if (li_.has_meta("JAC_" + name_ + "_SPARSITY_OUT", ind)) {
    return Sparsity::compressed(
      li_.meta_vector<casadi_int>("jac_" + name_ + "_SPARSITY_OUT", ind));
  } else {
    // Fall back to base class
    return FunctionInternal::get_jac_sparsity(oind, iind, symmetric);
  }
}

bool GenericExternal::get_diff_in(casadi_int i) {
  if (get_diff_in_) {
    // Query function exists
    return get_diff_in_(i);
  } else {
    // Fall back to base class
    return FunctionInternal::get_diff_in(i);
  }
}

bool GenericExternal::get_diff_out(casadi_int i) {
  if (get_diff_out_) {
    // Query function exists
    return get_diff_out_(i);
  } else {
    // Fall back to base class
    return FunctionInternal::get_diff_out(i);
  }
}

void External::init(const Dict& opts) {
  // Call the initialization method of the base class
  FunctionInternal::init(opts);

  casadi_assert(any_symbol_found(),
    "Could not find any function/symbol starting with '" + name_ + "_'. "
    "Make sure to read documentation of `external()` for proper usage.");

  // Reference counting?
  has_refcount_ = li_.has_function(name_ + "_incref");
  casadi_assert(has_refcount_==li_.has_function(name_ + "_decref"),
                        "External functions must provide functions for both increasing "
                        "and decreasing the reference count, or neither.");

  // Allocate work vectors
  casadi_int sz_arg=0, sz_res=0, sz_iw=0, sz_w=0;
  if (work_) {
    casadi_int flag = work_(&sz_arg, &sz_res, &sz_iw, &sz_w);
    casadi_assert(flag==0, "External: \"work\" failed");
  } else if (li_.has_meta(name_ + "_WORK")) {
    std::vector<casadi_int> v = li_.meta_vector<casadi_int>(name_ + "_WORK");
    casadi_assert_dev(v.size()==4);
    sz_arg = v[0];
    sz_res = v[1];
    sz_iw = v[2];
    sz_w = v[3];
  }

  // Work vectors
  alloc_arg(sz_arg);
  alloc_res(sz_res);
  alloc_iw(sz_iw);
  alloc_w(sz_w);
}

void GenericExternal::init(const Dict& opts) {
  // Call recursively
  External::init(opts);
}

void External::codegen_declarations(CodeGenerator& g) const {
  if (!li_.inlined(name_)) {
    g.add_external(signature(name_) + ";");
    if (checkout_) g.add_external("int " + name_ + "_checkout(void);");
    if (release_) g.add_external("void " + name_ + "_release(int mem);");
    if (incref_) g.add_external("void " + name_ + "_incref(void);");
    if (decref_) g.add_external("void " + name_ + "_decref(void);");
  }
}

void External::codegen_body(CodeGenerator& g) const {
  if (li_.inlined(name_)) {
    // Function body is inlined
    g << li_.body(name_) << "\n";
  } else {
    g << "if (" << name_ << "(arg, res, iw, w, mem)) return 1;\n";
  }
}

std::string External::codegen_mem_type() const {
  if (checkout_) return "nonempty";
  return "";
}

void External::codegen_checkout(CodeGenerator& g) const {
  if (checkout_) {
    g << "return " << name_ << "_checkout();\n";
  } else {
    g << "return 0;\n";
  }
}

void External::codegen_release(CodeGenerator& g) const {
  if (release_) {
    g << name_ << "_release(mem);\n";
  }
}

void External::codegen_incref(CodeGenerator& g) const {
  if (incref_) {
    g << name_ << "_incref();\n";
  }
}

void External::codegen_decref(CodeGenerator& g) const {
  if (decref_) {
    g << name_ << "_decref();\n";
  }
}

void External::codegen_init_mem(CodeGenerator& g) const {
  g << "return 0;\n";
}

void External::codegen_alloc_mem(CodeGenerator& g) const {
  g << "return 0;\n";
}

void External::codegen_free_mem(CodeGenerator& g) const {
}

bool External::has_jacobian() const {
  if (FunctionInternal::has_jacobian()) return true;
  return li_.has_function("jac_" + name_);
}

Function External
::get_jacobian(const std::string& name,
                  const std::vector<std::string>& inames,
                  const std::vector<std::string>& onames,
                  const Dict& opts) const {
  if (has_jacobian()) {
    return external(name, li_, opts);
  } else {
    return FunctionInternal::get_jacobian(name, inames, onames, opts);
  }
}

Function External
::get_forward(casadi_int nfwd, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames,
              const Dict& opts) const {
  // Consistency check
  casadi_int n=1;
  while (n<nfwd) n*=2;
  if (n!=nfwd || !has_forward(nfwd)) {
    // Inefficient code to be replaced later
    Function fwd1 = forward(1);
    return fwd1.map(name, "serial", nfwd, range(n_in_+n_out_), std::vector<casadi_int>(), opts);
    //casadi_error("Internal error: Refactoring needed, cf. #1055");
  }
  return external(name, li_, opts);
}

bool External::has_forward(casadi_int nfwd) const {
  return li_.has_function("fwd" + str(nfwd) + "_" + name_);
}

Function External
::get_reverse(casadi_int nadj, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames,
              const Dict& opts) const {
  // Consistency check
  casadi_int n=1;
  while (n<nadj) n*=2;
  if (n!=nadj || !has_reverse(nadj)) {
    // Inefficient code to be replaced later
    Function adj1 = reverse(1);
    return adj1.map(name, "serial", nadj, range(n_in_+n_out_), std::vector<casadi_int>(), opts);
    //casadi_error("Internal error: Refactoring needed, cf. #1055");
  }
  return external(name, li_, opts);
}

bool External::has_reverse(casadi_int nadj) const {
  return li_.has_function("adj" + str(nadj) + "_" + name_);
}

Function External::factory(const std::string& name,
                           const std::vector<std::string>& s_in,
                           const std::vector<std::string>& s_out,
                           const Function::AuxOut& aux,
                           const Dict& opts) const {
  // If not available, call base class function
  if (!li_.has_function(name)) {
    return FunctionInternal::factory(name, s_in, s_out, aux, opts);
  }

  // Retrieve function
  Function ret = external(name, li_, opts);

  // Replace colons in input names
  std::vector<std::string> s_io(s_in);
  for (std::string& s : s_io) replace(s.begin(), s.end(), ':', '_');

  // Inputs consistency checks
  casadi_assert(s_in.size() == ret.n_in(),
    "Inconsistent number of inputs. Expected " + str(s_in.size())+ "  "
    "(" + str(s_io) + "), got " + str(ret.n_in()) + ".");
  for (size_t i = 0; i < s_in.size(); ++i) {
    casadi_assert(s_io[i] == ret.name_in(i),
      "Inconsistent input name. Expected: " + str(s_io) + ", "
      "got " + ret.name_in(i) + " for input " + str(i));
  }

  // Replace colons in output names
  s_io = s_out;
  for (std::string& s : s_io) replace(s.begin(), s.end(), ':', '_');

  // Outputs consistency checks
  casadi_assert(s_out.size() == ret.n_out(),
    "Inconsistent number of outputs. Expected " + str(s_out.size()) + " "
    "(" + str(s_io) + "), got " + str(ret.n_out()) + ".");
  for (size_t i = 0; i < s_out.size(); ++i) {
    casadi_assert(s_io[i] == ret.name_out(i),
      "Inconsistent output name. Expected: " + str(s_io) + ", "
      "got " + ret.name_out(i) + " for output " + str(i));
  }

  return ret;
}

void External::serialize_body(SerializingStream &s) const {
  FunctionInternal::serialize_body(s);

  s.version("External", 1);
  s.pack("External::int_data", int_data_);
  s.pack("External::real_data", real_data_);
  s.pack("External::string_data", string_data_);
  s.pack("External::li", li_);
}

ProtoFunction* External::deserialize(DeserializingStream& s) {
  s.version("GenericExternal", 1);
  char type;
  s.unpack("GenericExternal::type", type);
  switch (type) {
    case 'g': return new GenericExternal(s);
    default:
      casadi_error("External::deserialize error");
  }
}

void GenericExternal::serialize_type(SerializingStream &s) const {
  FunctionInternal::serialize_type(s);
  s.version("GenericExternal", 1);
  s.pack("GenericExternal::type", 'g');
}

External::External(DeserializingStream & s) : FunctionInternal(s) {
  s.version("External", 1);
  s.unpack("External::int_data", int_data_);
  s.unpack("External::real_data", real_data_);
  s.unpack("External::string_data", string_data_);
  s.unpack("External::li", li_);

  External::init_external();
}

GenericExternal::GenericExternal(DeserializingStream& s) : External(s) {
  GenericExternal::init_external();
}

} // namespace casadi
