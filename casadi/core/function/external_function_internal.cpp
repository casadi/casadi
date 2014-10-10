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

ExternalFunctionInternal::ExternalFunctionInternal(const std::string& bin_name) :
    bin_name_(bin_name) {
#ifdef WITH_DL

  // Load the dll
#ifdef _WIN32
  handle_ = LoadLibrary(TEXT(bin_name_.c_str()));
  casadi_assert_message(handle_!=0, "ExternalFunctionInternal: Cannot open function: "
                        << bin_name_ << ". error code (WIN32): "<< GetLastError());

  initPtr init = (initPtr)GetProcAddress(handle_, TEXT("init"));
  if (init==0) throw CasadiException("ExternalFunctionInternal: no \"init\" found");
  getSparsityPtr getSparsity = (getSparsityPtr)GetProcAddress(handle_, TEXT("getSparsity"));
  if (getSparsity==0) throw CasadiException("ExternalFunctionInternal: no \"getSparsity\" found");
  evaluate_ = (evaluatePtr) GetProcAddress(handle_, TEXT("evaluateWrap"));
  if (evaluate_==0) throw CasadiException("ExternalFunctionInternal: no \"evaluateWrap\" found");

#else // _WIN32
  handle_ = dlopen(bin_name_.c_str(), RTLD_LAZY);
  casadi_assert_message(handle_!=0, "ExternalFunctionInternal: Cannot open function: "
                        << bin_name_ << ". error code: "<< dlerror());

  // reset error
  dlerror();

  // Load symbols
  initPtr init = (initPtr)dlsym(handle_, "init");
  if (dlerror()) throw CasadiException("ExternalFunctionInternal: no \"init\" found. "
                                       "Possible cause: If the function was generated from CasADi, "
                                       "make sure that it was compiled with a C compiler. If the "
                                       "function is C++, make sure to use extern \"C\" linkage.");
  getSparsityPtr getSparsity = (getSparsityPtr)dlsym(handle_, "getSparsity");
  if (dlerror()) throw CasadiException("ExternalFunctionInternal: no \"getSparsity\" found");
  evaluate_ = (evaluatePtr) dlsym(handle_, "evaluateWrap");
  if (dlerror()) throw CasadiException("ExternalFunctionInternal: no \"evaluateWrap\" found");
#endif // _WIN32

  // Initialize and get the number of inputs and outputs
  int n_in=-1, n_out=-1;
  int flag = init(&n_in, &n_out);
  if (flag) throw CasadiException("ExternalFunctionInternal: \"init\" failed");

  // Pass to casadi
  setNumInputs(n_in);
  setNumOutputs(n_out);

  // Get the sparsity pattern
  for (int i=0; i<n_in+n_out; ++i) {
    // Get sparsity from file
    int nrow, ncol, *colind, *row;
    flag = getSparsity(i, &nrow, &ncol, &colind, &row);
    if (flag) throw CasadiException("ExternalFunctionInternal: \"getSparsity\" failed");

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
      input(i) = Matrix<double>(sp, 0);
    } else {
      output(i-n_in) = Matrix<double>(sp, 0);
    }
  }

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

void ExternalFunctionInternal::evaluate() {
#ifdef WITH_DL
  int flag = evaluate_(getPtr(input_array_), getPtr(output_array_));
  if (flag) throw CasadiException("ExternalFunctionInternal: \"evaluate\" failed");
#endif // WITH_DL
}

void ExternalFunctionInternal::init() {
  // Call the init function of the base class
  FunctionInternal::init();

  // Get pointers to the inputs
  input_array_.resize(getNumInputs());
  for (int i=0; i<input_array_.size(); ++i)
    input_array_[i] = input(i).ptr();

  // Get pointers to the outputs
  output_array_.resize(getNumOutputs());
  for (int i=0; i<output_array_.size(); ++i)
    output_array_[i] = output(i).ptr();
}




} // namespace casadi

