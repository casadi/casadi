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
    // Default options
    string scheme = "forward1";
    double stepsize = 1e-8;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="scheme") {
        scheme = op.second.to_string();
      } else if (op.first=="stepsize") {
        stepsize = op.second;
      }
    }

    // Create instance
    Function ret;
    if (scheme=="forward1") {
      ret.assignNode(new Forward1(name, f, n, stepsize));
    } else if (scheme=="central1") {
      casadi_error("'central1' Not yet implemented");
    } else {
      casadi_error("No such scheme: '" + scheme + "'"
                   " Supported: 'forward1', 'central1'");
    }
    ret->construct(opts);
    return ret;
  }

  Derivative::Derivative(const std::string& name, const Function& f, int n, double h)
    : FunctionInternal(name), f_(f), n_(n), h_(h) {
  }

  Derivative::~Derivative() {
  }

  Options Derivative::options_
  = {{&FunctionInternal::options_},
     {{"stepsize",
       {OT_DOUBLE,
        "Perturbation size [default: 1e-8]"}},
      {"scheme",
       {OT_STRING,
        "Differencing scheme [default: 'forward1']"}}
     }
  };

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

  Sparsity Derivative::get_sparsity_in(int i) {
    int n_in = f_.n_in(), n_out = f_.n_out();
    if (i<n_in) {
      // Non-differentiated input
      return f_.sparsity_in(i);
    } else if (i<n_in+n_out) {
      // Non-differentiated output
      if (uses_output()) {
        return f_.sparsity_out(i-n_in);
      } else {
        return Sparsity(f_.size_out(i-n_in));
      }
    } else {
      // Seeds
      return repmat(f_.sparsity_in(i-n_in-n_out), 1, n_);
    }
  }

  Sparsity Derivative::get_sparsity_out(int i) {
    return repmat(f_.sparsity_out(i), 1, n_);
  }

  double Derivative::default_in(int ind) const {
    if (ind<f_.n_in()) {
      return f_.default_in(ind);
    } else {
      return 0;
    }
  }

  std::string Derivative::get_name_in(int i) {
    int n_in = f_.n_in(), n_out = f_.n_out();
    if (i<n_in) {
      return f_.name_in(i);
    } else if (i<n_in+n_out) {
      return "out_" + f_.name_out(i-n_in);
    } else {
      return "fwd_" + f_.name_in(i-n_in-n_out);
    }
  }

  std::string Derivative::get_name_out(int i) {
    return "fwd_" + f_.name_out(i);
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
      casadi_copy(f_res_pert, nnz, sens[j]);
      casadi_axpy(nnz, -1., f_res[j], sens[j]);
      casadi_scal(nnz, 1./h_, sens[j]);
      f_res_pert += nnz;
    }
  }

} // namespace casadi
