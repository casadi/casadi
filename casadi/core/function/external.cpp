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

  Function external(const string& name, const Dict& opts) {
    Function ret;
    ret.assignNode(External::create("./" + name + ".so", name));
    ret->construct(opts);
    return ret;
  }

  Function external(const string& name, const string& bin_name,
                    const Dict& opts) {
    Function ret;
    ret.assignNode(External::create(bin_name, name));
    ret->construct(opts);
    return ret;
  }

  Function external(const string& name, const Compiler& compiler,
                    const Dict& opts) {
    Function ret;
    ret.assignNode(External::create(compiler, name));
    ret->construct(opts);
    return ret;
  }

  LibInfo<std::string>::LibInfo(const std::string& bin_name)
    : bin_name_(bin_name), handle_(0) {
#ifdef WITH_DL
#ifdef _WIN32
    handle_ = LoadLibrary(TEXT(bin_name_.c_str()));
    casadi_assert_message(handle_!=0, "CommonExternal: Cannot open function: "
                          << bin_name_ << ". error code (WIN32): "<< GetLastError());
#else // _WIN32
    handle_ = dlopen(bin_name_.c_str(), RTLD_LAZY);
    casadi_assert_message(handle_!=0, "CommonExternal: Cannot open function: "
                          << bin_name_ << ". error code: "<< dlerror());
    // reset error
    dlerror();
#endif // _WIN32
#else // WITH_DL
    casadi_error("CommonExternal: WITH_DL  not activated");
#endif // WITH_DL
  }

  LibInfo<Compiler>::LibInfo(const Compiler& compiler)
    : compiler_(compiler) {
    if (meta().has("SYMBOLS")) meta_symbols_ = meta().to_set<std::string>("SYMBOLS");
  }

  void LibInfo<std::string>::clear() {
#ifdef WITH_DL
    // close the dll
#ifdef _WIN32
    if (handle_) FreeLibrary(handle_);
#else // _WIN32
    if (handle_) dlclose(handle_);
#endif // _WIN32
    handle_ = 0;
#endif // WITH_DL
  }

  template<typename FcnPtr>
  void LibInfo<std::string>::get(FcnPtr& fcnPtr, const std::string& sym) {
#ifdef WITH_DL
#ifdef _WIN32
    fcnPtr = (FcnPtr)GetProcAddress(handle_, TEXT(sym.c_str()));
#else // _WIN32
    fcnPtr = (FcnPtr)dlsym(handle_, sym.c_str());
    if (dlerror()) {
      fcnPtr=0;
      dlerror(); // Reset error flags
    }
#endif // _WIN32
#else // WITH_DL
#endif // WITH_DL
  }

  template<typename FcnPtr>
  void LibInfo<Compiler>::get(FcnPtr& fcnPtr, const std::string& sym) {
    fcnPtr = (FcnPtr)compiler_.get_function(sym);
  }

  bool LibInfo<Compiler>::has(const std::string& sym) {
    // Check if in meta information
    if (meta_symbols_.count(sym)) return true;

    // Convert to a dummy function pointer
    void (*fcn)(void);
    get(fcn, sym);
    return fcn!=0;
  }

  bool LibInfo<std::string>::has(const std::string& sym) {
    // Convert to a dummy function pointer
    void (*fcn)(void);
    get(fcn, sym);
    return fcn!=0;
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
  }

  External*
  External::create(const std::string& bin_name, const std::string& name) {
    return createGeneric(bin_name, name);
  }

  External*
  External::create(const Compiler& compiler, const std::string& name) {
    return createGeneric(compiler, name);
  }

  template<typename LibType>
  External* External::
  createGeneric(const LibType& libtype, const std::string& f_name) {
    // Structure with info about the library to be passed to the constructor
    LibInfo<LibType> li(libtype);

    if (li.has(f_name + "_simple")) {
      // Simplified, lower overhead external
      return new SimplifiedExternal<LibType>(f_name, li);
    } else {
      // Full information external
      return new GenericExternal<LibType>(f_name, li);
    }
  }

  template<typename LibType>
  CommonExternal<LibType>::CommonExternal(const std::string& name, const LibInfo<LibType>& li)
    : External(name), li_(li) {

    // Initialize and get the number of inputs and outputs
    // Function for creating an object instance
    string init_s = name + "_init";
    int n_int=0, n_real=0;
    li_.get(init_, init_s);
    if (init_!=0) {
      int flag = init_(0, 0, &n_int, &n_real);
      casadi_assert_message(flag==0, "CommonExternal: \"init\" failed");
    } else {
      // Fall back to reading meta information
      if (li_.meta().has(name + "_N_INT")) n_int=li_.meta().to_int(name + "_N_INT");
      if (li_.meta().has(name + "_N_REAL")) n_real=li_.meta().to_int(name + "_N_REAL");
    }
    int_data_.resize(n_int);
    real_data_.resize(n_real);
  }

  template<typename LibType>
  void CommonExternal<LibType>::init(const Dict& opts) {
    // Call recursively
    External::init(opts);

    // Get number of temporaries
    int sz_arg, sz_res, sz_iw, sz_w;
    work_t work;
    li_.get(work, name_ + "_work");
    if (work!=0) {
      int flag = work(&sz_arg, &sz_res, &sz_iw, &sz_w);
      casadi_assert_message(flag==0, "CommonExternal: \"work\" failed");
    } else {
      // No work vectors
      sz_arg = n_in();
      sz_res = n_out();
      sz_iw = sz_w = 0;
    }

    // Function for allocating memory, if any
    li_.get(allocmem_, name_ + "_alloc");

    // Function for freeing memory, if any
    li_.get(freemem_, name_ + "_free");

    alloc_arg(sz_arg);
    alloc_res(sz_res);
    alloc_iw(sz_iw);
    alloc_w(sz_w);
  }

  template<typename LibType>
  SimplifiedExternal<LibType>::SimplifiedExternal(const std::string& name,
                                                  const LibInfo<LibType>& li)
    : CommonExternal<LibType>(name, li) {
  }

  template<typename LibType>
  void SimplifiedExternal<LibType>::init(const Dict& opts) {
    // Call recursively
    CommonExternal<LibType>::init(opts);

    // Function for numerical evaluation
    string eval_s = this->name_ + "_simple";
    li_.get(eval_, eval_s);
    casadi_assert_message(eval_!=0, "SimplifiedExternal: no \"" + eval_s + "\" found");

    // Arrays for holding inputs and outputs
    this->alloc_w(this->n_in() + this->n_out());
  }

  template<typename LibType>
  size_t CommonExternal<LibType>::get_n_in() const {
    if (init_) {
      int n_in;
      int flag = init_(&n_in, 0, 0, 0);
      casadi_assert_message(flag==0, "CommonExternal: \"init\" failed");
      return n_in;
    } else if (li_.meta().has(name_ + "_N_IN")) {
      return li_.meta().to_int(name_ + "_N_IN");
    } else {
      // Fall back to base class
      return External::get_n_in();
    }
  }

  template<typename LibType>
  size_t CommonExternal<LibType>::get_n_out() const {
    if (init_) {
      int n_out;
      int flag = init_(0, &n_out, 0, 0);
      casadi_assert_message(flag==0, "CommonExternal: \"init\" failed");
      return n_out;
    } else if (li_.meta().has(name_ + "_N_OUT")) {
      return li_.meta().to_int(name_ + "_N_OUT");
    } else {
      // Fall back to base class
      return External::get_n_out();
    }
  }

  template<typename LibType>
  Sparsity CommonExternal<LibType>::get_sparsity_in(int ind) const {
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

  template<typename LibType>
  Sparsity CommonExternal<LibType>::get_sparsity_out(int ind) const {
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

  template<typename LibType>
  GenericExternal<LibType>::GenericExternal(const std::string& name, const LibInfo<LibType>& li)
    : CommonExternal<LibType>(name, li) {

    // Function for retrieving sparsities of inputs and outputs
    li_.get(this->sparsity_, name + "_sparsity");
  }

  External::External(const std::string& name)
    : FunctionInternal(name) {
  }

  External::~External() {
  }

  template<typename LibType>
  CommonExternal<LibType>::~CommonExternal() {
    if (freemem_) {
      freemem_(0);
    }
    li_.clear();
  }

  template<typename LibType>
  void GenericExternal<LibType>::init(const Dict& opts) {
    // Call recursively
    CommonExternal<LibType>::init(opts);

    // Function for numerical evaluation
    li_.get(eval_, this->name_);
  }

  template<typename LibType>
  void SimplifiedExternal<LibType>::
  simple(const double* arg, double* res) {
    casadi_assert_message(eval_!=0, "Numerical evaluation not possible");
    eval_(arg, res);
  }

  template<typename LibType>
  void GenericExternal<LibType>::
  eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
    casadi_assert_message(eval_!=0, "Numerical evaluation not possible");
    int flag = eval_(mem, arg, res, iw, w);
    if (flag) throw CasadiException("CommonExternal: \""+this->name_+"\" failed");
  }

  void External::addDependency(CodeGenerator& g) const {
    g.addExternal(signature(this->name_) + ";");
  }

  template<typename LibType>
  std::string SimplifiedExternal<LibType>::
  simple_call(const CodeGenerator& g,
              const std::string& arg, const std::string& res) const {
    // Create a function call
    stringstream ss;
    ss << this->name_ << "(" << arg << ", " << res << ")";
    return ss.str();
  }

  template<typename LibType>
  std::string GenericExternal<LibType>::
  generic_call(const CodeGenerator& g, const std::string& mem,
               const std::string& arg, const std::string& res,
               const std::string& iw, const std::string& w) const {

    // Create a function call
    stringstream ss;
    ss << this->name_ << "(" << mem << ", " << arg << ", " << res << ", "
       << iw << ", " << w << ")";
    return ss.str();
  }

  template<typename LibType>
  bool CommonExternal<LibType>::hasFullJacobian() const {
    if (FunctionInternal::hasFullJacobian()) return true;
    return const_cast<LibInfo<LibType>&>(li_).has(name_ + "_jac");
  }

  template<typename LibType>
  Function CommonExternal<LibType>
  ::getFullJacobian(const std::string& name, const Dict& opts) {
    if (hasFullJacobian()) {
      return external(name, li_, opts);
    } else {
      return FunctionInternal::getFullJacobian(name, opts);
    }
  }

  template<typename LibType>
  Function CommonExternal<LibType>
  ::get_forward(const std::string& name, int nfwd, Dict& opts) {
    // Consistency check
    int n=1;
    while (n<nfwd) n*=2;
    if (n!=nfwd || nfwd>get_n_forward()) {
      casadi_error("Internal error: Refactoring needed, cf. #1055");
    }
    return external(name, li_, opts);
  }

  template<typename LibType>
  int CommonExternal<LibType>::get_n_forward() const {
    // Will try 64, 32, 16, 8, 4, 2, 1 directions
    init_t fwd_init;
    for (int i=64; i>0; i/=2) {
      stringstream ss;
      ss << "fwd" << i << "_" << name_;
      const_cast<LibInfo<LibType>&>(li_).get(fwd_init, ss.str());
      if (fwd_init!=0) return i;
    }
    return 0;
  }

  template<typename LibType>
  Function CommonExternal<LibType>
  ::get_reverse(const std::string& name, int nadj, Dict& opts) {
    // Consistency check
    int n=1;
    while (n<nadj) n*=2;
    if (n!=nadj || nadj>get_n_reverse()) {
      casadi_error("Internal error: Refactoring needed, cf. #1055");
    }
    return external(name, li_, opts);
  }

  template<typename LibType>
  int CommonExternal<LibType>::get_n_reverse() const {
    // Will try 64, 32, 16, 8, 4, 2, 1 directions
    init_t adj_init;
    for (int i=64; i>0; i/=2) {
      stringstream ss;
      ss << "adj" << i << "_" << name_;
      const_cast<LibInfo<LibType>&>(li_).get(adj_init, ss.str());
      if (adj_init!=0) return i;
    }
    return 0;
  }

} // namespace casadi

