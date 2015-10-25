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
    n_in = n_out = sz_arg = sz_res = 0;
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
    n_in = n_out = sz_arg = sz_res = 0;
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

  ExternalFunctionInternal*
  ExternalFunctionInternal::create(const std::string& bin_name, const std::string& name) {
    return createGeneric(bin_name, name);
  }

  ExternalFunctionInternal*
  ExternalFunctionInternal::create(const Compiler& compiler, const std::string& name) {
    return createGeneric(compiler, name);
  }

  template<typename LibType>
  ExternalFunctionInternal* ExternalFunctionInternal::
  createGeneric(const LibType& libtype, const std::string& f_name) {
    // Structure with info about the library to be passed to the constructor
    LibInfo<LibType> li(libtype);

    // Function to retrieving number of inputs and outputs
    initPtr init;
    string init_s = f_name + "_init";
    li.get(init, init_s);
    if (init==0) {
      li.clear();
      casadi_error("CommonExternal: no \""+init_s+"\" found. "
                   "Possible cause: If the function was generated from CasADi, "
                   "make sure that it was compiled with a C compiler. If the "
                   "function is C++, make sure to use extern \"C\" linkage.");
    }

    // Initialize and get the number of inputs and outputs
    int f_type=0;
    int flag = init(&f_type, &li.n_in, &li.n_out, &li.sz_arg, &li.sz_res);
    if (flag) {
      li.clear();
      casadi_error("CommonExternal: \"init\" failed");
    }
    li.sz_arg = max(li.sz_arg, li.n_in);
    li.sz_res = max(li.sz_res, li.n_out);

    // Call the constructor
    switch (f_type) {
    case 0:
      // Simplified, lower overhead external
      return new SimplifiedExternal<LibType>(f_name, li);
    case 1:
      // Full information external
      return new GenericExternal<LibType>(f_name, li);
    default:
      li.clear();
      casadi_error("CommonExternal: Function type " << f_type << " not supported.");
    }
    return 0;
  }

  template<typename LibType>
  CommonExternal<LibType>::CommonExternal(const std::string& name, const LibInfo<LibType>& li)
    : ExternalFunctionInternal(name), li_(li) {
    ibuf_.resize(li.n_in);
    obuf_.resize(li.n_out);
    alloc_arg(li.sz_arg);
    alloc_res(li.sz_res);
  }

  template<typename LibType>
  SimplifiedExternal<LibType>::SimplifiedExternal(const std::string& name,
                                                  const LibInfo<LibType>& li)
    : CommonExternal<LibType>(name, li) {
    // Function for numerical evaluation
    li_.get(eval_, name);
    casadi_assert_message(eval_!=0, "SimplifiedExternal: no \""+name+"\" found");

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

  Sparsity ExternalFunctionInternal::get_sparsity_in(int ind) const {
    return get_sparsity(ind);
  }

  Sparsity ExternalFunctionInternal::get_sparsity_out(int ind) const {
    return get_sparsity(ind+n_in());
  }

  Sparsity ExternalFunctionInternal::get_sparsity(int ind) const {
    // Get sparsity from file
    int nrow, ncol;
    const int *colind, *row;
    casadi_assert(sparsity_!=0);
    int flag = sparsity_(ind, &nrow, &ncol, &colind, &row);
    casadi_assert_message(flag==0, "ExternalFunctionInternal: \"sparsity\" failed");

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

    // Get number of temporaries
    workPtr work;
    li_.get(work, name + "_work");
    if (work!=0) {
      int n_iw, n_w;
      int flag = work(&n_iw, &n_w);
      casadi_assert_message(flag==0, "CommonExternal: \"work\" failed");
      this->alloc_iw(n_iw);
      this->alloc_w(n_w);
    }
  }

  ExternalFunctionInternal::ExternalFunctionInternal(const std::string& name)
    : FunctionInternal(name) {
  }

  ExternalFunctionInternal::~ExternalFunctionInternal() {
  }

  template<typename LibType>
  CommonExternal<LibType>::~CommonExternal() {
    li_.clear();
  }

  template<typename LibType>
  void SimplifiedExternal<LibType>::evalD(const double** arg, double** res,
                                          int* iw, double* w) {
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
  void GenericExternal<LibType>::evalD(const double** arg, double** res,
                                       int* iw, double* w) {
    int flag = eval_(arg, res, iw, w);
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
  std::string SimplifiedExternal<LibType>::generateCall(const CodeGenerator& g,
                                                        const std::string& arg,
                                                        const std::string& res) const {
    // Create a function call
    stringstream ss;
    ss << this->name_ << "(" << arg << ", " << res << ")";
    return ss.str();
  }

  template<typename LibType>
  std::string GenericExternal<LibType>::
  generateCall(const CodeGenerator& g,
               const std::string& arg, const std::string& res,
               const std::string& iw, const std::string& w) const {

    // Create a function call
    stringstream ss;
    ss << this->name_ << "(" << arg << ", " << res << ", " << iw << ", " << w << ")";
    return ss.str();
  }

  template<typename LibType>
  int GenericExternal<LibType>::scalarSparsity(int i, int *n_row, int *n_col,
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
    initPtr jac_init;
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
  ::getDerForward(const std::string& name, int nfwd, Dict& opts) {
    // Consistency check
    int n=1;
    while (n<nfwd) n*=2;
    if (n!=nfwd || nfwd>numDerForward()) {
      casadi_error("Internal error: Refactoring needed, cf. #1055");
    }
    return Function::external(name, li_, opts);
  }

  template<typename LibType>
  int CommonExternal<LibType>::numDerForward() const {
    // Will try 64, 32, 16, 8, 4, 2, 1 directions
    initPtr fwd_init;
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
  ::getDerReverse(const std::string& name, int nadj, Dict& opts) {
    // Consistency check
    int n=1;
    while (n<nadj) n*=2;
    if (n!=nadj || nadj>numDerReverse()) {
      casadi_error("Internal error: Refactoring needed, cf. #1055");
    }
    return Function::external(name, li_, opts);
  }

  template<typename LibType>
  int CommonExternal<LibType>::numDerReverse() const {
    // Will try 64, 32, 16, 8, 4, 2, 1 directions
    initPtr adj_init;
    for (int i=64; i>0; i/=2) {
      stringstream ss;
      ss << "adj" << i << "_" << name_;
      const_cast<LibInfo<LibType>&>(li_).get(adj_init, ss.str());
      if (adj_init!=0) return i;
    }
    return 0;
  }

} // namespace casadi

