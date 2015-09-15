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

#include "callback.hpp"
#include "../functor.hpp"
#include "mx_function.hpp"

using namespace std;

namespace casadi {

std::vector<DMatrix> Callback2::operator()(const std::vector<DMatrix>& arg) {
  DMatrix out = arg[0];
  std::vector<DMatrix> ret;
  ret.push_back(out*2);
  return ret;
}

Callback2::Callback2() {

}


Callback2::~Callback2() {

}

Function finiteDiffGenerator(Function& fcn, int ndir, void* user_data) {
  // The step size to take
  DMatrix eps = * reinterpret_cast<double *>(user_data);

  // Obtain the symbols for nominal inputs/outputs
  std::vector<MX> nominal_in  = fcn.symbolicInput();
  std::vector<MX> nominal_out = fcn.symbolicOutput();

  // A growing list of inputs to the returned derivative function
  std::vector<MX> der_ins = nominal_in;
  der_ins.insert(der_ins.end(),  nominal_out.begin(),  nominal_out.end());

  // A growing list of outputs to the returned derivative function
  std::vector<MX> der_outs;

  for (int k=0;k<ndir;++k) {
     std::vector<MX> seeds  = fcn.symbolicInput();
     std::vector<MX> perturbed_in(seeds.size());
     for (int i=0;i<perturbed_in.size();++i) {
       perturbed_in[i] = nominal_in[i] + eps*seeds[i];
     }
     std::vector<MX> perturbed_out  = fcn(perturbed_in);
     std::vector<MX> sens(perturbed_out.size());
     for (int i=0;i<perturbed_out.size();++i) {
       sens[i] = (perturbed_out[i]-nominal_out[i])/eps;
     }
     der_ins.insert(der_ins.end(),   seeds.begin(),  seeds.end());
     der_outs.insert(der_outs.end(), sens.begin(),   sens.end());
  }
  MXFunction ret("finite_diff", der_ins, der_outs);
  return ret;
}

Function Callback2::create() {
  Function ret;
  CallbackFunctionInternal* node = new CallbackFunctionInternal(*this);
  ret.assignNode(node);
  ret.setOption("name", name());
  ret.setOption("custom_forward", DerivativeGenerator(finiteDiffGenerator) );

  ret.setOption(options());
  node->eps_ = ret.getOption("fin_diff_eps");

  ret.setOption("user_data", reinterpret_cast<void *>(&node->eps_));
  ret.init();
  return ret;
}

Function DerivativeGenerator2::operator()(Function& fcn, int ndir) {
  return fcn->getDerForward("der", ndir, Dict());
  casadi_error("This virtual method must be implemented");
}

Function DerivativeGenerator2::original(Function& fcn, int ndir, bool fwd) {
  if (fwd) {
    return fcn->getDerForward("der", ndir, Dict());
  } else {
    return fcn->getDerReverse("der", ndir, Dict());
  }
}

DerivativeGenerator DerivativeGenerator2::create() {
  DerivativeGenerator ret;
  ret.assignNode(new DerivativeGeneratorInternal2(*this));
  return ret;
}

DerivativeGenerator2::DerivativeGenerator2() {
}

DerivativeGenerator2::~DerivativeGenerator2() {

}

DerivativeGeneratorInternal2::DerivativeGeneratorInternal2(
    DerivativeGenerator2 &callback) : callback_(callback) {

}

DerivativeGeneratorInternal2::~DerivativeGeneratorInternal2() {

}


CallbackFunctionInternal::CallbackFunctionInternal(
    Callback2 &callback) : callback_(callback) {

  addOption("fin_diff_eps", OT_REAL, 1e-7, "eps used for finite differences");

  ibuf_.resize(callback_.nIn());
  obuf_.resize(callback_.nOut());

  for (int k=0;k<ibuf_.size();k++) {
    input(k) = DMatrix::zeros(callback_.inputSparsity(k));
  }

  for (int k=0;k<obuf_.size();k++) {
    output(k) = DMatrix::zeros(callback_.outputSparsity(k));
  }

}

CallbackFunctionInternal::~CallbackFunctionInternal() {

}


void CallbackFunctionInternal::evalD(const double** arg,
                               double** res, int* iw, double* w) {
    // Number of inputs and outputs
    int num_in = nIn();
    int num_out = nOut();

    std::vector<DMatrix> inputs(num_in);

    // Pass the inputs to the function
    for (int i=0; i<num_in; ++i) {
      inputs[i] = DMatrix::zeros(input(i).sparsity());
      if (arg[i] != 0) {
        inputs[i].setNZ(arg[i]);
      } else {
        inputs[i].set(0.);
      }
    }

    std::vector<DMatrix> outputs = callback_(inputs);

    // Get the outputs
    for (int i=0; i<num_out; ++i) {
      if (res[i] != 0) outputs[i].getNZ(res[i]);
    }
  }

void CallbackFunctionInternal::init() {
  FunctionInternal::init();
}



} // namespace casadi

