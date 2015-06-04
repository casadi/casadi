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


#include "external_function_internal.hpp"
#include "../std_vector_tools.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

namespace casadi {
  using namespace std;

   ExternalFunctionInternal::handle_t
   ExternalFunctionInternal::getHandle(const std::string& bin_name) {
     handle_t handle;
#ifdef WITH_DL
#ifdef _WIN32
    handle = LoadLibrary(TEXT(bin_name.c_str()));
    casadi_assert_message(handle!=0, "ExternalFunctionInternal: Cannot open function: "
                          << bin_name << ". error code (WIN32): "<< GetLastError());
#else // _WIN32
    handle = dlopen(bin_name.c_str(), RTLD_LAZY);
    casadi_assert_message(handle!=0, "ExternalFunctionInternal: Cannot open function: "
                          << bin_name << ". error code: "<< dlerror());
    // reset error
    dlerror();
#endif // _WIN32
#else // WITH_DL
    casadi_error("ExternalFunctionInternal: WITH_DL  not activated");
#endif // WITH_DL
    return handle;
   }

  template<typename FcnPtr>
  void ExternalFunctionInternal::getSym(FcnPtr& fcnPtr, handle_t handle,
                                        const std::string& sym) {
#ifdef WITH_DL
#ifdef _WIN32
    fcnPtr = (FcnPtr)GetProcAddress(handle, TEXT(sym.c_str()));
#else // _WIN32
    fcnPtr = (FcnPtr)dlsym(handle, sym.c_str());
    if (dlerror()) {
      fcnPtr=0;
      dlerror(); // Reset error flags
    }
#endif // _WIN32
#else // WITH_DL
#endif // WITH_DL
  }

  void ExternalFunctionInternal::freeHandle(handle_t& handle) {
#ifdef WITH_DL
    // close the dll
#ifdef _WIN32
    if (handle) {
      FreeLibrary(handle);
      handle=0;
    }
#else // _WIN32
    if (handle) {
      dlclose(handle);
      handle=0;
    }
#endif // _WIN32
#endif // WITH_DL
  }

  ExternalFunctionInternal*
  ExternalFunctionInternal::create(const std::string& bin_name, const std::string& f_name) {
    // Structure with info about the library to be passed to the constructor
    LibInfo li = {bin_name, f_name, 0, 0, 0, 0, 0};

    // Load the dll and get the handle
    li.handle = getHandle(bin_name);

    // Function to retrieving number of inputs and outputs
    initPtr init;
    string init_s = f_name + "_init";
    getSym(init, li.handle, init_s);
    if (init==0) {
      freeHandle(li.handle);
      casadi_error("ExternalFunctionInternal: no \""+init_s+"\" found. "
                   "Possible cause: If the function was generated from CasADi, "
                   "make sure that it was compiled with a C compiler. If the "
                   "function is C++, make sure to use extern \"C\" linkage.");
    }

    // Initialize and get the number of inputs and outputs
    int f_type=0;
    int flag = init(&f_type, &li.n_in, &li.n_out, &li.sz_arg, &li.sz_res);
    if (flag) {
      freeHandle(li.handle);
      casadi_error("ExternalFunctionInternal: \"init\" failed");
    }
    li.sz_arg = max(li.sz_arg, li.n_in);
    li.sz_res = max(li.sz_res, li.n_out);

    // Call the constructor
    switch (f_type) {
    case 0:
      // Simplified, lower overhead external
      return new SimplifiedExternal(li);
    case 1:
      // Full information external
      return new GenericExternal(li);
    default:
      freeHandle(li.handle);
      casadi_error("ExternalFunctionInternal: Function type " << f_type << " not supported.");
    }
    return 0;
  }

  ExternalFunctionInternal::ExternalFunctionInternal(const LibInfo& li)
    : bin_name_(li.bin_name), f_name_(li.f_name), handle_(li.handle) {
    ibuf_.resize(li.n_in);
    obuf_.resize(li.n_out);
    alloc_arg(li.sz_arg);
    alloc_res(li.sz_res);

    // Set name TODO(@jaeandersson): remove 'name' from options
    setOption("name", f_name_);
  }

  SimplifiedExternal::SimplifiedExternal(const LibInfo& li) : ExternalFunctionInternal(li) {
    // Function for numerical evaluation
    getSym(eval_, handle_, f_name_);
    casadi_assert_message(eval_!=0, "SimplifiedExternal: no \""+f_name_+"\" found");

    // All inputs are scalars
    for (int i=0; i<nIn(); ++i) {
      input(i) = 0;
    }

    // All outputs are scalars
    for (int i=0; i<nOut(); ++i) {
      output(i) = 0;
    }

    // Arrays for holding inputs and outputs
    alloc_w(nIn() + nOut());
  }

  GenericExternal::GenericExternal(const LibInfo& li) : ExternalFunctionInternal(li) {
    // Function for numerical evaluation
    getSym(eval_, handle_, f_name_);
    casadi_assert_message(eval_!=0, "GenericExternal: no \""+f_name_+"\" found");

    // Function for retrieving sparsities of inputs and outputs
    sparsityPtr sparsity;
    getSym(sparsity, handle_, f_name_ + "_sparsity");
    if (sparsity==0) {
      // Fall back to scalar sparsity
      sparsity = scalarSparsity;
    }

    // Get the sparsity patterns
    for (int i=0; i<nIn()+nOut(); ++i) {
      // Get sparsity from file
      int nrow, ncol;
      const int *colind, *row;
      int flag = sparsity(i, &nrow, &ncol, &colind, &row);
      casadi_assert_message(flag==0, "ExternalFunctionInternal: \"sparsity\" failed");

      // Col offsets
      vector<int> colindv(colind, colind+ncol+1);

      // Number of nonzeros
      int nnz = colindv.back();

      // Rows
      vector<int> rowv(row, row+nnz);

      // Sparsity
      Sparsity sp(nrow, ncol, colindv, rowv);

      // Save to inputs/outputs
      if (i<nIn()) {
        input(i) = Matrix<double>::zeros(sp);
      } else {
        output(i-nIn()) = Matrix<double>::zeros(sp);
      }
    }

    // Get number of temporaries
    workPtr work;
    getSym(work, handle_, f_name_ + "_work");
    if (work!=0) {
      int n_iw, n_w;
      int flag = work(&n_iw, &n_w);
      casadi_assert_message(flag==0, "ExternalFunctionInternal: \"work\" failed");
      alloc_iw(n_iw);
      alloc_w(n_w);
    }
  }

  ExternalFunctionInternal* ExternalFunctionInternal::clone() const {
    throw CasadiException("Error ExternalFunctionInternal cannot be cloned");
  }

  ExternalFunctionInternal::~ExternalFunctionInternal() {
    freeHandle(handle_);
  }

  void SimplifiedExternal::evalD(const double** arg, double** res,
                                 int* iw, double* w) {
    // Copy arguments to input buffers
    const double* arg1=w;
    for (int i=0; i<nIn(); ++i) {
      *w++ = arg[i] ? *arg[i] : 0;
    }

    // Evaluate
    eval_(arg1, w);

    // Get outputs
    for (int i=0; i<nOut(); ++i) {
      if (res[i]) *res[i] = *w;
      ++w;
    }
  }

  void GenericExternal::evalD(const double** arg, double** res,
                                       int* iw, double* w) {
    int flag = eval_(arg, res, iw, w);
    if (flag) throw CasadiException("ExternalFunctionInternal: \""+f_name_+"\" failed");
  }

  void ExternalFunctionInternal::init() {
    // Call the init function of the base class
    FunctionInternal::init();
  }

  void SimplifiedExternal::addDependency(CodeGenerator& g) const {
    g.addExternal("void " + f_name_ + "(const real_t* arg, real_t* res);");
  }

  void GenericExternal::addDependency(CodeGenerator& g) const {
    g.addExternal("int " + f_name_ + "(const real_t** arg, real_t** res, int* iw, real_t* w);");
  }

  std::string SimplifiedExternal::generateCall(const CodeGenerator& g,
                                               const std::string& arg,
                                               const std::string& res) const {
    // Create a function call
    stringstream ss;
    ss << f_name_ << "(" << arg << ", " << res << ")";
    return ss.str();
  }

  std::string GenericExternal::generateCall(const CodeGenerator& g,
                                            const std::string& arg, const std::string& res,
                                            const std::string& iw, const std::string& w) const {
    // Create a function call
    stringstream ss;
    ss << f_name_ << "(" << arg << ", " << res << ", " << iw << ", " << w << ")";
    return ss.str();
  }

  int GenericExternal::scalarSparsity(int i, int *n_row, int *n_col,
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
} // namespace casadi

