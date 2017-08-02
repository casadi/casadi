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

  Function CentralDiff::create(const std::string& name, int n, const Dict& opts) {
    return Function::create(new CentralDiff(name, n), opts);
  }

  CentralDiff::CentralDiff(const std::string& name, int n)
    : FunctionInternal(name), n_(n) {
  }

  CentralDiff::~CentralDiff() {
  }

  Options CentralDiff::options_
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

  void CentralDiff::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Default options
    h_ = 1e-8;
    h2_ = 1e-3;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="stepsize") {
        h_ = op.second;
      } else if (op.first=="second_order_stepsize") {
        h2_ = op.second;
      } else if (op.first=="scheme") {
        casadi_warning("Option 'scheme' currently ignored");
      }
    }

    // Allocate work vector for (perturbed) inputs and outputs
    n_x_ = derivative_of_.nnz_in();
    n_f_ = derivative_of_.nnz_out();
    alloc_w(3 * n_x_, true); // x, x0, v
    alloc_w(3 * n_f_, true); // f, f0, Jv

    // Allocate sufficient temporary memory for function evaluation
    alloc(derivative_of_);
  }

  Sparsity CentralDiff::get_sparsity_in(int i) {
    int n_in = derivative_of_.n_in(), n_out = derivative_of_.n_out();
    if (i<n_in) {
      // Non-differentiated input
      return derivative_of_.sparsity_in(i);
    } else if (i<n_in+n_out) {
      // Non-differentiated output
      return derivative_of_.sparsity_out(i-n_in);
    } else {
      // Seeds
      return repmat(derivative_of_.sparsity_in(i-n_in-n_out), 1, n_);
    }
  }

  Sparsity CentralDiff::get_sparsity_out(int i) {
    return repmat(derivative_of_.sparsity_out(i), 1, n_);
  }

  double CentralDiff::default_in(int ind) const {
    if (ind<derivative_of_.n_in()) {
      return derivative_of_.default_in(ind);
    } else {
      return 0;
    }
  }

  size_t CentralDiff::get_n_in() {
    return derivative_of_.n_in() + derivative_of_.n_out() + derivative_of_.n_in();
  }

  size_t CentralDiff::get_n_out() {
    return derivative_of_.n_out();
  }

  std::string CentralDiff::get_name_in(int i) {
    int n_in = derivative_of_.n_in(), n_out = derivative_of_.n_out();
    if (i<n_in) {
      return derivative_of_.name_in(i);
    } else if (i<n_in+n_out) {
      return "out_" + derivative_of_.name_out(i-n_in);
    } else {
      return "fwd_" + derivative_of_.name_in(i-n_in-n_out);
    }
  }

  std::string CentralDiff::get_name_out(int i) {
    return "fwd_" + derivative_of_.name_out(i);
  }

  Function CentralDiff::get_forward(int nfwd, const std::string& name,
                                   const std::vector<std::string>& inames,
                                   const std::vector<std::string>& onames,
                                   const Dict& opts) const {
    Dict opts_mod = opts;
    opts_mod["stepsize"] = h2_;
    return Function::create(new CentralDiff(name, nfwd), opts_mod);
  }

  void CentralDiff::eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
    // Memory structure
    Mem m_tmp;
    Mem *m = &m_tmp;
    m->n_x = n_x_;
    m->n_f = n_f_;
    m->h = h_;

    // Shorthands
    int n_in = derivative_of_.n_in(), n_out = derivative_of_.n_out();

    // Non-differentiated input
    m->x = w;
    for (int j=0; j<n_in; ++j) {
      const int nnz = derivative_of_.nnz_in(j);
      casadi_copy(*arg++, nnz, w);
      w += nnz;
    }

    // Non-differentiated output
    m->f = w;
    for (int j=0; j<n_out; ++j) {
      const int nnz = derivative_of_.nnz_out(j);
      casadi_copy(*arg++, nnz, w);
      w += nnz;
    }

    // Forward seeds
    const double** seed = arg; arg += n_in;

    // Forward sensitivities
    double** sens = res; res += n_out;

    // Vector with which Jacobian is multiplied
    m->v = w; w += m->n_x;
    casadi_fill(m->v, m->n_x, 0.);

    // Jacobian-times-vector product
    m->Jv = w; w += m->n_f;

    // Unperturbed input
    m->x0 = w; w += m->n_x;

    // Unperturbed output
    m->f0 = w; w += m->n_f;

    // Setup arg, res for calling f
    double* x1 = m->x;
    for (int j=0; j<n_in; ++j) {
      arg[j] = x1;
      x1 += derivative_of_.nnz_in(j);
    }
    double* f1 = m->f;
    for (int j=0; j<n_out; ++j) {
      res[j] = f1;
      f1 += derivative_of_.nnz_out(j);
    }


    // For each derivative direction
    for (int i=0; i<n_; ++i) {
      // Copy seeds to v
      double* v1 = m->v;
      for (int j=0; j<n_in; ++j) {
        int nnz = derivative_of_.nnz_in(j);
        if (seed[j]) casadi_copy(seed[j] + nnz*i, nnz, v1);
        v1 += nnz;
      }

      // Reset counter
      m->n_calls = 0;

      // Call reverse communication algorithm
      while (central_differences(m)) {
        derivative_of_(arg, res, iw, w, 0);
      }

      // Gather sensitivities
      double* Jv1 = m->Jv;
      for (int j=0; j<n_out; ++j) {
        int nnz = derivative_of_.nnz_out(j);
        if (sens[j]) casadi_copy(Jv1, nnz, sens[j] + i*nnz);
        Jv1 += nnz;
      }
    }
  }

  bool CentralDiff::central_differences(Mem* m) {
    bool ret = true;
    switch (m->n_calls) {
      case 0:
      // Backup x and f
      casadi_copy(m->x, m->n_x, m->x0);
      casadi_copy(m->f, m->n_f, m->f0);

      // Perturb x, positive direction
      casadi_axpy(m->n_x, m->h/2, m->v, m->x);
      break; // first function call
      case 1:

      // Save result, perturb in negative direction
      casadi_copy(m->f, m->n_f, m->Jv);
      casadi_copy(m->x0, m->n_x, m->x);
      casadi_axpy(m->n_x, -m->h/2, m->v, m->x);
      break; // second function call
      case 2:

      // Calculate finite difference approximation
      casadi_axpy(m->n_f, -1., m->f, m->Jv);
      casadi_scal(m->n_f, 1/m->h, m->Jv);

      // Restore x and f
      casadi_copy(m->x0, m->n_x, m->x);
      casadi_copy(m->f0, m->n_f, m->f);
      ret = false;
    }

    // Increase function call counter
    if (ret) m->n_calls++;
    return ret;
  }

} // namespace casadi
