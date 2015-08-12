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


#include "mapaccum_internal.hpp"
#include "mx_function.hpp"

using namespace std;

namespace casadi {

  MapAccumInternal::MapAccumInternal(const Function& f, int n)
    : f_(f), n_(n) {

    casadi_assert(f.nIn()>=1);
    casadi_assert(f.nOut()>=1);

    casadi_assert_message(f.inputSparsity(0)==f.outputSparsity(0),
                            "First input and output must have matching sparsity." <<
                            "Got " << f.inputSparsity(0).dimString() << " and " <<
                            f.outputSparsity(0).dimString() << ".");

    // Give a name
    setOption("name", "unnamed_mapaccum");
  }

  MapAccumInternal::~MapAccumInternal() {
  }

  void MapAccumInternal::init() {

    int num_in = f_.nIn(), num_out = f_.nOut();

    // Initialize the functions, get input and output sparsities
    // Input and output sparsities

    ibuf_.resize(num_in);
    obuf_.resize(num_out);

    for (int i=0;i<num_in;++i) {
      // Sparsity of the original input
      Sparsity in_sp = f_.input(i).sparsity();

      if (i>=1) in_sp = repmat(in_sp, 1, n_);

      // Allocate space for input
      input(i) = DMatrix::zeros(in_sp);
    }

    for (int i=0;i<num_out;++i) {
      // Sparsity of the original output
      Sparsity out_sp = f_.output(i).sparsity();

      if (i>=1) out_sp = repmat(out_sp, 1, n_);

      // Allocate space for output
      output(i) = DMatrix::zeros(out_sp);
    }

    step_in_.resize(num_in, 0);
    step_out_.resize(num_out, 0);

    for (int i=0;i<num_in;++i) {
       step_in_[i] = f_.input(i).nnz();
    }

    for (int i=0;i<num_out;++i) {
       step_out_[i] = f_.output(i).nnz();
    }


    // Call the initialization method of the base class
    FunctionInternal::init();

    alloc_w(f_.sz_w()+step_in_[0]);
    alloc_iw(f_.sz_iw());
    alloc_arg(2*f_.sz_arg());
    alloc_res(2*f_.sz_res());
  }

  template<typename T, typename R>
  void MapAccumInternal::evalGen(const T** arg, T** res, int* iw, T* w,
    void (FunctionInternal::*eval)(const T** arg, T** res, int* iw, T* w),
    R reduction) {
    int num_in = f_.nIn(), num_out = f_.nOut();

    const T** arg1 = arg+f_.sz_arg();
    T** res1 = res+f_.sz_res();

    // Reset the accumulator
    T* accum = w+f_.sz_w();
    std::copy(arg[0], arg[0]+step_in_[0], accum);
    for (int i=0; i<n_; ++i) {

      // Set the function inputs
      arg1[0] = accum;
      for (int j=1; j<num_in; ++j) {
        arg1[j] = (arg[j]==0) ? 0: arg[j]+i*step_in_[j];
      }
      // Set the function outputs
      res1[0] = res[0];
      for (int j=1; j<num_out; ++j) {
        res1[j] = (res[j]==0)? 0: res[j]+i*step_out_[j];
      }
      // Evaluate the function
      FunctionInternal *f = f_.operator->();
      (f->*eval)(arg1, res1, iw, w);

      // Copy the accumulator
      std::copy(res1[0], res1[0]+step_in_[0], accum);

    }
  }

  void MapAccumInternal::evalD(const double** arg, double** res,
                                int* iw, double* w) {
    evalGen<double>(arg, res, iw, w, &FunctionInternal::evalD, std::plus<double>());
  }

  void MapAccumInternal::evalSX(const SXElement** arg, SXElement** res,
                                int* iw, SXElement* w) {
  }

  void MapAccumInternal::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
  }

  Function MapAccumInternal
  ::getDerForward(const std::string& name, int nfwd, const Dict& opts) {

  }

  Function MapAccumInternal
  ::getDerReverse(const std::string& name, int nadj, const Dict& opts) {

  }

  void MapAccumInternal::generateDeclarations(CodeGenerator& g) const {
    f_->addDependency(g);
  }

  void MapAccumInternal::generateBody(CodeGenerator& g) const {

  }

  inline string name(const Function& f) {
    if (f.isNull()) {
      return "NULL";
    } else {
      return f.getOption("name");
    }
  }

  void MapAccumInternal::print(ostream &stream) const {
    stream << "MapAccum(" << name(f_) << ", " << n_ << ")";
  }

} // namespace casadi
