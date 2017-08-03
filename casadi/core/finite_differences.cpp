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


#include "finite_differences.hpp"

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
    n_z_ = derivative_of_.nnz_in();
    n_f_ = derivative_of_.nnz_out();
    alloc_w(2 * n_, true); // m->x, m->x0
    alloc_w(2 * n_f_, true); // m->f, m->f0
    alloc_w(n_*n_f_, true); // m->J
    alloc_w(n_z_, true); // z

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
    // Shorthands
    int n_in = derivative_of_.n_in(), n_out = derivative_of_.n_out();

    // Non-differentiated input
    const double** z0 = arg; arg += n_in;

    // Non-differentiated output
    double* f = w;
    for (int j=0; j<n_out; ++j) {
      const int nnz = derivative_of_.nnz_out(j);
      casadi_copy(*arg++, nnz, w);
      w += nnz;
    }

    // Forward seeds
    const double** seed = arg; arg += n_in;

    // Forward sensitivities
    double** sens = res; res += n_out;

    // Memory structure
    casadi_central_diff_mem<double> m_tmp, *m = &m_tmp;
    m->n_x = n_;
    m->n_f = n_f_;
    m->h = h_;
    m->f = f;
    m->x = w; w += n_;
    m->J = w; w += n_f_*n_;
    m->x0 = w; w += n_;
    m->f0 = w; w += n_f_;

    // central_diff only sees a function with n_ inputs, initialized at 0
    casadi_fill(m->x, n_, 0.);

    // Setup arg, res for calling function
    double* z = w;
    for (int j=0; j<n_in; ++j) {
      arg[j] = w;
      w += derivative_of_.nnz_in(j);
    }
    for (int j=0; j<n_out; ++j) {
      res[j] = f;
      f += derivative_of_.nnz_out(j);
    }

    // Call reverse communication algorithm
    m->next = 0;
    while (casadi_central_diff(m)) {
      double* z1 = z;
      for (int j=0; j<n_in; ++j) {
        int nnz = derivative_of_.nnz_in(j);
        casadi_copy(z0[j], nnz, z1);
        casadi_mv_dense(seed[j], nnz, n_, m->x, z1, false);
        z1 += nnz;
      }
      derivative_of_(arg, res, iw, w, 0);
    }

    // Gather sensitivities
    for (int i=0; i<n_; ++i) {
      for (int j=0; j<n_out; ++j) {
        int nnz = derivative_of_.nnz_out(j);
        if (sens[j]) casadi_copy(m->J, nnz, sens[j] + i*nnz);
        m->J += nnz;
      }
    }
  }

  void CentralDiff::codegen_declarations(CodeGenerator& g) const {
    derivative_of_->add_dependency(g);
  }

  void CentralDiff::codegen_body(CodeGenerator& g) const {
    // Shorthands
    int n_in = derivative_of_.n_in(), n_out = derivative_of_.n_out();

    g.comment("Non-differentiated input");
    g.local("z0", "const real_t", "**");
    g << "z0 = arg; arg += " << n_in << ";\n";

    g.comment("Non-differentiated output");
    g.local("f", "real_t", "*");
    g << "f = w;\n";
    for (int j=0; j<n_out; ++j) {
      const int nnz = derivative_of_.nnz_out(j);
      g << g.copy("*arg++", nnz, "w") << " w += " << nnz << ";\n";
    }

    g.comment("Forward seeds");
    g.local("seed", "const real_t", "**");
    g << "seed = arg; arg += " << n_in << ";\n";

    g.comment("Forward sensitivities");
    g.local("sens", "real_t", "**");
    g << "sens = res; res += " << n_out << ";\n";

    g.comment("Memory structure");
    g.local("m_tmp", "central_diff_mem");
    g.local("m", "central_diff_mem", "*");
    g << "m = &m_tmp;\n"
      << "m->n_x = " << n_ << ";\n"
      << "m->n_f = " << n_f_ << ";\n"
      << "m->h = " << h_ << ";\n"
      << "m->f = f;\n"
      << "m->x = w; w += " << n_ << ";\n"
      << "m->J = w; w += " << n_f_*n_ << ";\n"
      << "m->x0 = w; w += " << n_ << ";\n"
      << "m->f0 = w; w += " << n_f_ << ";\n";
    g << g.fill("m->x", n_, "0.") << "\n";

    g.comment("Setup arg, res for calling function");
    g.local("z", "real_t", "*");
    g << "z = w;\n";
    for (int j=0; j<n_in; ++j) {
      g << "arg[" << j << "] = w; w += " << derivative_of_.nnz_in(j) << ";\n";
    }
    for (int j=0; j<n_out; ++j) {
      g << "res[" << j << "] = f; f += " << derivative_of_.nnz_out(j) << ";\n";
    }

    g.comment("Invoke reverse communication algorithm");
    g << "m->next = 0;\n";
    g << "while (" << g.central_diff("m") << ") {\n";
    g.local("z1", "real_t", "*");
    g << "z1 = z;\n";
    for (int j=0; j<n_in; ++j) {
      int nnz = derivative_of_.nnz_in(j);
      g << g.copy("z0[" + to_string(j) + "]", nnz, "z1") << "\n"
        << g.mv("seed[" + to_string(j) + "]", nnz, n_, "m->x", "z1", false) << "\n"
        << "z1 += " << nnz << ";\n";
    }
    g << "if (" << g(derivative_of_, "arg", "res", "iw", "w") << ") return 1;\n"
      << "}\n"; // while (...)

    g.comment("Gather sensitivities");
    g.local("i", "int");
    g << "for (i=0; i<" << n_ << "; ++i) {" << "\n";
    for (int j=0; j<n_out; ++j) {
      int nnz = derivative_of_.nnz_out(j);
      string s = "sens[" + to_string(j) + "]";
      g << "if (" << s << ") " << g.copy("m->J", nnz, s + "+i*" + to_string(nnz)) << "\n"
        << "m->J += " << nnz << ";\n";
    }
    g << "}\n"; // for (i=0, ...)
  }


} // namespace casadi
