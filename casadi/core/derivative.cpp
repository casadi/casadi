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


#include "derivative.hpp"

using namespace std;

namespace casadi {

  Function Derivative::create(const std::string& name, const Function& f,
                              int n, const Dict& opts) {
    Function ret;
    ret.assignNode(new Forward1(name, f, n));
    ret->construct(opts);
    return ret;
  }

  Derivative::Derivative(const std::string& name, const Function& f, int n)
    : FunctionInternal(name), f_(f), n_(n) {
  }

  Derivative::~Derivative() {
  }

  void Derivative::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Allocate work vector for perturbed inputs
    alloc_w(n_calls() * f_.nnz_in(), true);
    alloc_w(n_calls() * f_.nnz_out(), true);

    // Work vectors for seeds/sensitivities
    alloc_arg(f_.n_in(), true);
    alloc_res(f_.n_out(), true);

    // Allocate sufficient temporary memory for function evaluation
    alloc(f_);
  }

  void Derivative::eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
    // Shorthands
    int n_in = f_.n_in(), n_out = f_.n_out(), n_calls = this->n_calls();

    // Non-differentiated input
    const double** f_arg = arg; arg += n_in;

    // Non-differentiated output
    const double** f_res = arg; arg += n_out;

    // Forward seeds
    const double** seed = arg; arg += n_in;

    // Forward sensitivities
    double** sens = res; res += n_out;

    // Copy sensitivitity arguments to temporary vectors to allow modification
    copy_n(seed, n_in, arg);
    seed = arg; arg += n_in;
    copy_n(sens, n_out, res);
    sens = res; res += n_out;

    // Work vectors for perturbed inputs and outputs
    double* f_arg_pert = w; w += n_calls * f_.nnz_in();
    double* f_res_pert = w; w += n_calls * f_.nnz_out();

    // For each derivative direction
    for (int i=0; i<n_; ++i) {
      // Perturb function argument (depends on differentiation algorithm)
      perturb(f_arg, f_arg_pert, seed);

      // Function evaluation
      double* f_arg_pert1 = f_arg_pert;
      double* f_res_pert1 = f_res_pert;
      for (int c=0; c<n_calls; ++c) {
        // Function inputs
        for (int j=0; j<n_in; ++j) {
          arg[j] = f_arg_pert1;
          f_arg_pert1 += f_.nnz_in(j);
        }
        // Function outputs
        for (int j=0; j<n_out; ++j) {
          res[j] = f_res_pert1;
          f_res_pert1 += f_.nnz_out(j);
        }
        // Call function
        f_(arg, res, iw, w, 0);
      }

      // Calculate finite difference approximation
      finalize(f_res, f_res_pert, sens);

      // Proceed to the next direction
      for (int j=0; j<n_in; ++j) if (seed[j]) seed[j] += f_.nnz_in(j);
      for (int j=0; j<n_out; ++j) if (sens[j]) sens[j] += f_.nnz_out(j);
    }
  }

  void Forward1::perturb(const double** f_arg, double* f_arg_pert, const double** seed) const {
    int n_in = f_.n_in();
    for (int j=0; j<n_in; ++j) {
      const int nnz = f_.nnz_in(j);
      casadi_copy(f_arg[j], nnz, f_arg_pert);
      casadi_axpy(nnz, h_, seed[j], f_arg_pert);
      f_arg_pert += nnz;
    }
  }

  void Forward1::finalize(const double** f_res, const double* f_res_pert, double** sens) const {
    int n_out = f_.n_out();
    for (int j=0; j<n_out; ++j) {
      const int nnz = f_.nnz_out(j);
      if (sens[j]) {
        casadi_copy(f_res_pert, nnz, sens[j]);
        casadi_axpy(nnz, -1., f_res[j], sens[j]);
        casadi_scal(nnz, 1./h_, sens[j]);
      }
      f_res_pert += nnz;
    }
  }

} // namespace casadi
