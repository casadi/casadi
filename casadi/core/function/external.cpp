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


#include "external.hpp"
#include "../std_vector_tools.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

namespace casadi {
  using namespace std;

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
    fcnPtr = (FcnPtr)compiler_.getFunction(sym);
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

    // Can the function be evaluated using the simple API?
    simple_t simple;
    li.get(simple, f_name + "_simple");
    if (simple!=0) {
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

    // Function for creating an object instance
    init_t init;
    string init_s = name + "_init";
    li_.get(init, init_s);
    if (init==0) {
      li_.clear();
      casadi_error("CommonExternal: no \"" + init_s + "\" found. "
                   "Possible cause: If the function was generated from CasADi, "
                   "make sure that it was compiled with a C compiler. If the "
                   "function is C++, make sure to use extern \"C\" linkage.");
    }

    // Initialize and get the number of inputs and outputs
    int flag = init(&mem_, &n_in_, &n_out_, 0);
    casadi_assert_message(flag==0, "CommonExternal: \"init\" failed");

    // Get number of temporaries
    int sz_arg, sz_res, sz_iw, sz_w;
    work_t work;
    li_.get(work, name + "_work");
    if (work!=0) {
      flag = work(mem_, &sz_arg, &sz_res, &sz_iw, &sz_w);
      casadi_assert_message(flag==0, "CommonExternal: \"work\" failed");
    } else {
      // No work vectors
      sz_arg = n_in_;
      sz_res = n_out_;
      sz_iw = sz_w = 0;
    }

    // Function for freeing memory, if any
    li_.get(freemem_, name + "_free");

    ibuf_.resize(n_in_);
    obuf_.resize(n_out_);
    alloc_arg(sz_arg);
    alloc_res(sz_res);
    alloc_iw(sz_iw);
    alloc_w(sz_w);
  }

  template<typename LibType>
  SimplifiedExternal<LibType>::SimplifiedExternal(const std::string& name,
                                                  const LibInfo<LibType>& li)
    : CommonExternal<LibType>(name, li) {
    // Function for numerical evaluation
    string eval_s = name + "_simple";
    li_.get(eval_, eval_s);
    casadi_assert_message(eval_!=0, "SimplifiedExternal: no \"" + eval_s + "\" found");

    // All inputs are scalars
    for (int i=0; i<this->n_in(); ++i) {
      this->input(i) = 0;
    }

    // All outputs are scalars
    for (int i=0; i<this->n_out(); ++i) {
      this->output(i) = 0;
    }

    // Arrays for holding inputs and outputs
    this->alloc_w(this->n_in() + this->n_out());
  }

  Sparsity External::get_sparsity_in(int ind) const {
    return get_sparsity(ind);
  }

  Sparsity External::get_sparsity_out(int ind) const {
    return get_sparsity(ind+n_in());
  }

  template<typename LibType>
  Sparsity CommonExternal<LibType>::get_sparsity(int ind) const {
    // Get sparsity from file
    int nrow, ncol;
    const int *colind, *row;
    casadi_assert(sparsity_!=0);
    int flag = sparsity_(mem_, ind, &nrow, &ncol, &colind, &row);
    casadi_assert_message(flag==0, "External: \"sparsity\" failed");

    // Col offsets
    vector<int> colindv(colind, colind+ncol+1);

    // Number of nonzeros
    int nnz = colindv.back();

    // Rows
    vector<int> rowv(row, row+nnz);

    // Sparsity
    return Sparsity(nrow, ncol, colindv, rowv);
  }

  template<typename LibType>
  GenericExternal<LibType>::GenericExternal(const std::string& name, const LibInfo<LibType>& li)
    : CommonExternal<LibType>(name, li) {
    // Function for numerical evaluation
    li_.get(eval_, name);
    casadi_assert_message(eval_!=0, "GenericExternal: no \""+name+"\" found");

    // Function for retrieving sparsities of inputs and outputs
    li_.get(this->sparsity_, name + "_sparsity");
    if (this->sparsity_==0) {
      // Fall back to scalar sparsity
      this->sparsity_ = scalarSparsity;
    }

    // Get input sparsities
    for (int i=0; i<this->n_in(); ++i) {
      this->input(i) = Matrix<double>::zeros(this->get_sparsity_in(i));
    }

    // Get output sparsities
    for (int i=0; i<this->n_out(); ++i) {
      this->output(i) = Matrix<double>::zeros(this->get_sparsity_out(i));
    }
  }

  External::External(const std::string& name)
    : FunctionInternal(name) {
  }

  External::~External() {
  }

  template<typename LibType>
  CommonExternal<LibType>::~CommonExternal() {
    if (freemem_) {
      freemem_(mem_);
    }
    li_.clear();
  }

  template<typename LibType>
  void SimplifiedExternal<LibType>::
  eval(const double** arg, double** res, int* iw, double* w, void* mem) {
    // Copy arguments to input buffers
    const double* arg1=w;
    for (int i=0; i<this->n_in(); ++i) {
      *w++ = arg[i] ? *arg[i] : 0;
    }

    // Evaluate
    eval_(arg1, w);

    // Get outputs
    for (int i=0; i<this->n_out(); ++i) {
      if (res[i]) *res[i] = *w;
      ++w;
    }
  }

  template<typename LibType>
  void GenericExternal<LibType>::
  eval(const double** arg, double** res, int* iw, double* w, void* mem) {
    int flag = eval_(arg, res, iw, w, mem);
    if (flag) throw CasadiException("CommonExternal: \""+this->name_+"\" failed");
  }

  template<typename LibType>
  void SimplifiedExternal<LibType>::addDependency(CodeGenerator& g) const {
    g.addExternal("void " + this->name_ + "(const real_t* arg, real_t* res);");
  }

  template<typename LibType>
  void GenericExternal<LibType>::addDependency(CodeGenerator& g) const {
    g.addExternal("int " + this->name_ + "(const real_t** arg, real_t** res, int* iw, real_t* w);");
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
  int GenericExternal<LibType>::scalarSparsity(void* mem, int i, int *n_row, int *n_col,
                                               const int **colind, const int **row) {
    // Dense scalar
    static const int colind_scalar[] = {0, 1};
    static const int row_scalar[] = {0};

    // Return to user
    if (n_row) *n_row=1;
    if (n_col) *n_col=1;
    if (colind) *colind=colind_scalar;
    if (row) *row=row_scalar;
    return 0;
  }

  template<typename LibType>
  bool CommonExternal<LibType>::hasFullJacobian() const {
    if (FunctionInternal::hasFullJacobian()) return true;
    // Init function for Jacobian?
    init_t jac_init;
    const_cast<LibInfo<LibType>&>(li_).get(jac_init, name_ + "_jac_init");
    return jac_init!=0;
  }

  template<typename LibType>
  Function CommonExternal<LibType>
  ::getFullJacobian(const std::string& name, const Dict& opts) {
    if (hasFullJacobian()) {
      return Function::external(name, li_, opts);
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
    return Function::external(name, li_, opts);
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
    return Function::external(name, li_, opts);
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

