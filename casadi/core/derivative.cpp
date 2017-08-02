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

  Function Derivative::create(const std::string& name, const Dict& opts) {
    // Default options
    string scheme = "central";
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
    if (scheme=="forward") {
      return Function::create(new Forward(name, stepsize), opts);
    } else if (scheme=="central") {
      return Function::create(new Central(name, stepsize), opts);
    } else {
      casadi_error("No such scheme: '" + scheme + "'"
                   " Supported: 'central', 'forward'");
    }
  }

  Derivative::Derivative(const std::string& name, double h)
    : FunctionInternal(name), h_(h) {
  }

  Derivative::~Derivative() {
  }

  Options Derivative::options_
  = {{&FunctionInternal::options_},
     {{"stepsize",
       {OT_DOUBLE,
        "Perturbation size [default: 1e-8]"}},
      {"second_order_stepsize",
       {OT_DOUBLE,
        "Second order perturbation size [default: 1e-3]"}},
      {"scheme",
       {OT_STRING,
        "Differencing scheme [default: 'central']"}}
     }
  };

  void Derivative::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Default options
    h2_ = 1e-3;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="second_order_stepsize") {
        h2_ = op.second;
      }
    }

    // Allocate work vector for perturbed inputs
    alloc_w(n_calls() * f().nnz_in(), true);
    alloc_w(n_calls() * f().nnz_out(), true);

    // Allocate sufficient temporary memory for function evaluation
    alloc(f());
  }

  Sparsity Derivative::get_sparsity_in(int i) {
    int n_in = derivative_of_.n_in();
    if (i<n_in) {
      // Non-differentiated input
      return derivative_of_.sparsity_in(i);
    } else {
      // Non-differentiated output
      if (uses_output()) {
        return derivative_of_.sparsity_out(i-n_in);
      } else {
        return Sparsity(derivative_of_.size_out(i-n_in));
      }
    }
  }

  Sparsity Derivative::get_sparsity_out(int i) {
    return Sparsity::dense(derivative_of_.numel_out(), derivative_of_.numel_in());
  }

  double Derivative::default_in(int ind) const {
    if (ind<derivative_of_.n_in()) {
      return derivative_of_.default_in(ind);
    } else {
      return 0;
    }
  }

  size_t Derivative::get_n_in() {
    return derivative_of_.n_in() + derivative_of_.n_out();
  }

  size_t Derivative::get_n_out() {
    return 1;
  }

  std::string Derivative::get_name_in(int i) {
    int n_in = derivative_of_.n_in();
    if (i<n_in) {
      return derivative_of_.name_in(i);
    } else {
      return "out_" + derivative_of_.name_out(i-n_in);
    }
  }

  std::string Derivative::get_name_out(int i) {
    return "jac";
  }

  void Derivative::eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
    // Shorthands
    int n_in = derivative_of_.n_in(), n_out = derivative_of_.n_out(), n_calls = this->n_calls();

    // Non-differentiated input
    const double** f_arg = arg; arg += n_in;

    // Non-differentiated output
    const double** f_res = arg; arg += n_out;

    // Jacobian
    double* jac = res[0];

    // Work vectors for perturbed inputs and outputs
    double* f_arg_pert = w; w += n_calls * f().nnz_in();
    double* f_res_pert = w; w += n_calls * f().nnz_out();

    // For each derivative direction
    int numel_in = derivative_of_.numel_in();
    int numel_out = derivative_of_.numel_out();
    for (int i=0; i<numel_in; ++i) {
      // Perturb function argument (depends on differentiation algorithm)
      perturb(i, f_arg, f_arg_pert);

      // Function evaluation
      double* f_arg_pert1 = f_arg_pert;
      double* f_res_pert1 = f_res_pert;
      for (int c=0; c<n_calls; ++c) {
        // Function inputs
        for (int j=0; j<n_in; ++j) {
          arg[j] = f_arg_pert1;
          f_arg_pert1 += f().nnz_in(j);
        }
        // Function outputs
        for (int j=0; j<n_out; ++j) {
          res[j] = f_res_pert1;
          f_res_pert1 += f().nnz_out(j);
        }
        // Call function
        f()(arg, res, iw, w, 0);
      }

      // Calculate finite difference approximation
      finalize(i, f_res, f_res_pert, jac);
      jac += numel_out;
    }
  }

  void Forward::perturb(int i, const double** f_arg, double* f_arg_pert) const {
    int n_in = derivative_of_.n_in();
    for (int sign=0; sign<2; ++sign) {
      for (int j=0; j<n_in; ++j) {
        const int nnz = derivative_of_.nnz_in(j);
        casadi_copy(f_arg[j], nnz, f_arg_pert);
        if (sign) {
          casadi_axpy(nnz, h_, seed[j], f_arg_pert);
        }
        f_arg_pert += nnz;
      }
    }
  }

  void Forward::finalize(const double** f_res, const double* f_res_pert, double* jac) const {
    const double* f_res_pert1 = f_res_pert + derivative_of_.nnz_out();
    int n_out = derivative_of_.n_out();
    for (int j=0; j<n_out; ++j) {
      const int nnz = derivative_of_.nnz_out(j);
      casadi_copy(f_res_pert, nnz, sens[j]);
      f_res_pert += nnz;
      casadi_axpy(nnz, -1., f_res_pert1, sens[j]);
      casadi_scal(nnz, -1/h_, sens[j]);
      f_res_pert1 += nnz;
    }
  }

  void Central::perturb(const double** f_arg, double* f_arg_pert, const double** seed) const {
    int n_in = derivative_of_.n_in();
    for (int sign=0; sign<2; ++sign) {
      for (int j=0; j<n_in; ++j) {
        const int nnz = derivative_of_.nnz_in(j);
        casadi_copy(f_arg[j], nnz, f_arg_pert);
        casadi_axpy(nnz, sign ? -h_/2 : h_/2, seed[j], f_arg_pert);
        f_arg_pert += nnz;
      }
    }
  }

  void Central::finalize(const double** f_res, const double* f_res_pert, double* jac) const {
    const double* f_res_pert1 = f_res_pert + derivative_of_.nnz_out();
    int n_out = derivative_of_.n_out();
    for (int j=0; j<n_out; ++j) {
      const int nnz = derivative_of_.nnz_out(j);
      casadi_copy(f_res_pert, nnz, sens[j]);
      f_res_pert += nnz;
      casadi_axpy(nnz, -1., f_res_pert1, sens[j]);
      casadi_scal(nnz, 1/h_, sens[j]);
      f_res_pert1 += nnz;
    }
  }

} // namespace casadi
