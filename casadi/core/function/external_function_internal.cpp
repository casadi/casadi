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

  ExternalFunctionInternal::ExternalFunctionInternal(const string& bin_name, const string& f_name) :
    bin_name_(bin_name), f_name_(f_name) {
#ifdef WITH_DL

    // Names of the functions we want to access
    string init_s = f_name_ + "_init";
    string sparsity_s = f_name_ + "_sparsity";
    string work_s = f_name_ + "_work";
    bool has_work = true;

    // Load the dll
#ifdef _WIN32
    handle_ = LoadLibrary(TEXT(bin_name_.c_str()));
    casadi_assert_message(handle_!=0, "ExternalFunctionInternal: Cannot open function: "
                          << bin_name_ << ". error code (WIN32): "<< GetLastError());

    // Function to retrieving number of inputs and outputs
    initPtr init = (initPtr)GetProcAddress(handle_, TEXT(init_s.c_str()));
    if (init==0) throw CasadiException("ExternalFunctionInternal: no \""+init_s+"\" found");

    // Function for numerical evaluation
    eval_ = (evalPtr) GetProcAddress(handle_, TEXT(f_name_.c_str()));
    if (eval_==0) throw CasadiException("ExternalFunctionInternal: no \""+f_name_+"\" found");

    // Function for retrieving sparsities of inputs and outputs
    sparsityPtr sparsity = (sparsityPtr)GetProcAddress(handle_, TEXT(sparsity_s.c_str()));
    if (sparsity==0) {
      // Fall back to scalar sparsity
      sparsity = scalarSparsity;
    }

    // Function for retriving sizes of required work vectors
    workPtr work = (workPtr)GetProcAddress(handle_, TEXT(work_s.c_str()));
    if (work==0) {
      // No work vectors needed
      has_work = false;
    }
#else // _WIN32
    handle_ = dlopen(bin_name_.c_str(), RTLD_LAZY);
    casadi_assert_message(handle_!=0, "ExternalFunctionInternal: Cannot open function: "
                          << bin_name_ << ". error code: "<< dlerror());

    // reset error
    dlerror();

    // Function to retrieving number of inputs and outputs
    initPtr init = (initPtr)dlsym(handle_, init_s.c_str());
    casadi_assert_message(!dlerror(), "ExternalFunctionInternal: no \""+init_s+"\" found. "
                          "Possible cause: If the function was generated from CasADi, "
                          "make sure that it was compiled with a C compiler. If the "
                          "function is C++, make sure to use extern \"C\" linkage.");

    // Function for numerical evaluation
    eval_ = (evalPtr) dlsym(handle_, f_name_.c_str());
    if (dlerror()) throw CasadiException("ExternalFunctionInternal: no \""+f_name_+"\" found");

    // Function for retrieving sparsities of inputs and outputs
    sparsityPtr sparsity = (sparsityPtr)dlsym(handle_, sparsity_s.c_str());
    if (dlerror()) {
      // Fall back to scalar sparsity
      sparsity = scalarSparsity;
      // Reset error flags
      dlerror();
    }

    // Function for retriving sizes of required work vectors
    workPtr work = (workPtr)dlsym(handle_, work_s.c_str());
    if (dlerror()) {
      // No work vectors needed
      has_work = false;
      // Reset error flags
      dlerror();
    }

#endif // _WIN32
    // Initialize and get the number of inputs and outputs
    int f_type=-1, n_in=-1, n_out=-1, n_arg=-1, n_res=-1;
    int flag = init(&f_type, &n_in, &n_out, &n_arg, &n_res);
    if (flag) throw CasadiException("ExternalFunctionInternal: \"init\" failed");
    casadi_assert(f_type>=0 && n_in>=0 && n_out>=0 && n_arg>=0);
    casadi_assert(n_in<=n_arg && n_out<=n_res);

    // Pass to casadi
    ibuf_.resize(n_in);
    obuf_.resize(n_out);
    alloc_arg(n_arg);
    alloc_res(n_res);

    // Get the sparsity pattern
    for (int i=0; i<n_in+n_out; ++i) {
      // Get sparsity from file
      int nrow, ncol;
      const int *colind, *row;
      flag = sparsity(i, &nrow, &ncol, &colind, &row);
      if (flag) throw CasadiException("ExternalFunctionInternal: \"sparsity\" failed");

      // Col offsets
      vector<int> colindv(colind, colind+ncol+1);

      // Number of nonzeros
      int nnz = colindv.back();

      // Rows
      vector<int> rowv(row, row+nnz);

      // Sparsity
      Sparsity sp = Sparsity(nrow, ncol, colindv, rowv);

      // Save to inputs/outputs
      if (i<n_in) {
        input(i) = Matrix<double>::zeros(sp);
      } else {
        output(i-n_in) = Matrix<double>::zeros(sp);
      }
    }

    // Get number of temporaries
    int n_iw, n_w;
    if (has_work) {
      flag = work(&n_iw, &n_w);
      if (flag) throw CasadiException("ExternalFunctionInternal: \"work\" failed");
    } else {
      n_iw=n_w=0;
    }
    alloc_iw(n_iw);
    alloc_w(n_w);
#else // WITH_DL
    throw CasadiException("WITH_DL  not activated");
#endif // WITH_DL

    // Default name
    setOption("name", f_name);
  }

  ExternalFunctionInternal* ExternalFunctionInternal::clone() const {
    throw CasadiException("Error ExternalFunctionInternal cannot be cloned");
  }

  ExternalFunctionInternal::~ExternalFunctionInternal() {
#ifdef WITH_DL
    // close the dll
#ifdef _WIN32
    if (handle_) FreeLibrary(handle_);
#else // _WIN32
    if (handle_) dlclose(handle_);
#endif // _WIN32
#endif // WITH_DL
  }

  void ExternalFunctionInternal::evalD(const double** arg, double** res,
                                       int* iw, double* w) {
#ifdef WITH_DL
    int flag = eval_(arg, res, iw, w);
    if (flag) throw CasadiException("ExternalFunctionInternal: \""+f_name_+"\" failed");
#endif // WITH_DL
  }

  void ExternalFunctionInternal::init() {
    // Call the init function of the base class
    FunctionInternal::init();
  }

  void ExternalFunctionInternal::generateDeclarations(CodeGenerator& g) const {
    // Declare function (definition in separate file)
    g.body
      << "/* Defined in " << bin_name_ << " */" << endl
      << "int " << f_name_ << "(const real_t* const* arg, real_t* const* res, "
      << "int* iw, real_t* w);" << endl << endl;
  }

  void ExternalFunctionInternal::generateBody(CodeGenerator& g) const {
    g.body
      << "  int flag = " << f_name_ << "(arg, res, iw, w);" << endl
      << "  if (flag) return flag;" << endl;
  }

  int ExternalFunctionInternal::scalarSparsity(int i, int *n_row, int *n_col,
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

