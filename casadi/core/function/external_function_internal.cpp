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
    string narg_s = f_name_ + "_narg";
    string sparsity_s = f_name_ + "_sparsity";
    string work_s = f_name_ + "_work";

    // Load the dll
#ifdef _WIN32
    handle_ = LoadLibrary(TEXT(bin_name_.c_str()));
    casadi_assert_message(handle_!=0, "ExternalFunctionInternal: Cannot open function: "
                          << bin_name_ << ". error code (WIN32): "<< GetLastError());

    nargPtr narg = (nargPtr)GetProcAddress(handle_, TEXT(narg_s.c_str()));
    if (narg==0) throw CasadiException("ExternalFunctionInternal: no \""+narg_s+"\" found");
    sparsityPtr sparsity = (sparsityPtr)GetProcAddress(handle_, TEXT(sparsity_s.c_str()));
    if (sparsity==0)
      throw CasadiException("ExternalFunctionInternal: no \""+sparsity_s+"\" found");
    eval_ = (evalPtr) GetProcAddress(handle_, TEXT(f_name_.c_str()));
    if (eval_==0) throw CasadiException("ExternalFunctionInternal: no \""+f_name_+"\" found");
    workPtr work = (workPtr)GetProcAddress(handle_, TEXT(work_s.c_str()));
    if (work==0) throw CasadiException("ExternalFunctionInternal: no \""+work_s+"\" found");

#else // _WIN32
    handle_ = dlopen(bin_name_.c_str(), RTLD_LAZY);
    casadi_assert_message(handle_!=0, "ExternalFunctionInternal: Cannot open function: "
                          << bin_name_ << ". error code: "<< dlerror());

    // reset error
    dlerror();

    // Load symbols
    nargPtr narg = (nargPtr)dlsym(handle_, narg_s.c_str());
    casadi_assert_message(!dlerror(), "ExternalFunctionInternal: no \""+narg_s+"\" found. "
                          "Possible cause: If the function was generated from CasADi, "
                          "make sure that it was compiled with a C compiler. If the "
                          "function is C++, make sure to use extern \"C\" linkage.");
    sparsityPtr sparsity = (sparsityPtr)dlsym(handle_, sparsity_s.c_str());
    if (dlerror()) throw CasadiException("ExternalFunctionInternal: no \""+sparsity_s+"\" found");
    eval_ = (evalPtr) dlsym(handle_, f_name_.c_str());
    if (dlerror()) throw CasadiException("ExternalFunctionInternal: no \""+f_name_+"\" found");
    workPtr work = (workPtr)dlsym(handle_, work_s.c_str());
    if (dlerror()) throw CasadiException("ExternalFunctionInternal: no \""+work_s+"\" found");

#endif // _WIN32
    // Initialize and get the number of inputs and outputs
    int n_in=-1, n_out=-1;
    int flag = narg(&n_in, &n_out);
    if (flag) throw CasadiException("ExternalFunctionInternal: \"narg\" failed");

    // Pass to casadi
    setNumInputs(n_in);
    setNumOutputs(n_out);

    // Get the sparsity pattern
    for (int i=0; i<n_in+n_out; ++i) {
      // Get sparsity from file
      int nrow, ncol, *colind, *row;
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
    int ni, nr;
    flag = work(&ni, &nr);
    if (flag) throw CasadiException("ExternalFunctionInternal: \"work\" failed");
    ni_ = static_cast<size_t>(ni);
    nr_ = static_cast<size_t>(nr);
#else // WITH_DL
    throw CasadiException("WITH_DL  not activated");
#endif // WITH_DL

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

  void ExternalFunctionInternal::evalD(cp_double* arg, p_double* res,
                                       int* itmp, double* rtmp) {
#ifdef WITH_DL
    int flag = eval_(arg, res, itmp, rtmp);
    if (flag) throw CasadiException("ExternalFunctionInternal: \""+f_name_+"\" failed");
#endif // WITH_DL
  }

  void ExternalFunctionInternal::init() {
    // Call the init function of the base class
    FunctionInternal::init();
  }

  void ExternalFunctionInternal::nTmp(size_t& ni, size_t& nr) const {
    ni=ni_;
    nr=nr_;
  }

  void ExternalFunctionInternal::generateDeclarations(CodeGenerator& gen) const {
    // Declare function (definition in separate file)
    gen.functions
      << "/* Defined in " << bin_name_ << " */" << endl
      << "int " << f_name_ << "(const " << gen.real_t << "* const* arg, " << gen.real_t
      << "* const* res, int* iii, " << gen.real_t << "* w);" << endl << endl;
  }

  void ExternalFunctionInternal::generateBody(CodeGenerator& gen) const {
    gen.functions
      << "  int flag = " << f_name_ << "(arg, res, iii, w);" << endl
      << "  if (flag) return flag;" << endl;
  }

} // namespace casadi

