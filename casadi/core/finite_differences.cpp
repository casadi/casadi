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
    alloc_w(3 * n_, true); // m->x, m->x0, v
    alloc_w(3 * n_f_, true); // m->f, m->f0, m->Jv
    alloc_w(3 * n_z_, true); // z, z0, vz

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
    const double* z0 = w;
    for (int j=0; j<n_in; ++j) {
      const int nnz = derivative_of_.nnz_in(j);
      casadi_copy(*arg++, nnz, w);
      w += nnz;
    }

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
    m->n_x = 1;
    m->n_f = n_f_;
    m->h = h_;
    m->f = f;

    // Input
    m->x = w; w += n_;
    casadi_fill(m->x, n_, 0.);

    // Vector with which Jacobian is multiplied
    double* vz = w; w += n_z_;
    m->v = w; w += n_;
    casadi_fill(m->v, n_, 1.);

    // Jacobian-times-vector product
    m->Jv = w; w += m->n_f;

    // Unperturbed input
    m->x0 = w; w += n_;

    // Unperturbed output
    m->f0 = w; w += m->n_f;

    // Setup arg, res for calling function
    double* z = w;
    for (int j=0; j<n_in; ++j) {
      arg[j] = w;
      w += derivative_of_.nnz_in(j);
    }
    double* f1 = f;
    for (int j=0; j<n_out; ++j) {
      res[j] = f1;
      f1 += derivative_of_.nnz_out(j);
    }

    // For each derivative direction
    for (int i=0; i<n_; ++i) {
      // Copy seeds to v
      double* v1 = vz;
      for (int j=0; j<n_in; ++j) {
        int nnz = derivative_of_.nnz_in(j);
        if (seed[j]) casadi_copy(seed[j] + nnz*i, nnz, v1);
        v1 += nnz;
      }

      // Reset counter
      m->n_calls = 0;

      // Call reverse communication algorithm
      while (casadi_central_diff(m)) {
        casadi_copy(z0, n_z_, z);
        casadi_axpy(n_z_, m->x[0], vz, z);
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

  void CentralDiff::codegen_declarations(CodeGenerator& g) const {
    derivative_of_->add_dependency(g);
  }

  void CentralDiff::codegen_body(CodeGenerator& g) const {
    // Shorthands
    int n_in = derivative_of_.n_in(), n_out = derivative_of_.n_out();

    g.comment("Memory structure");
    g.local("m_tmp", "central_diff_mem");
    g.local("m", "central_diff_mem", "*");
    g << "m = &m_tmp;\n"
      << "m->n_x = " << n_z_ << ";\n"
      << "m->n_f = " << n_f_ << ";\n"
      << "m->h = " << h_ << ";\n";

    g.comment("Non-differentiated input");
    g << "m->x = w;\n";
    for (int j=0; j<n_in; ++j) {
      const int nnz = derivative_of_.nnz_in(j);
      g << g.copy("*arg++", nnz, "w") << " w += " << nnz << ";\n";
    }

    g.comment("Non-differentiated output");
    g << "m->f = w;\n";
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

    g.comment("Vector with which Jacobian is multiplied");
    g << "m->v = w; w += " << n_z_ << ";\n";
    g << g.fill("m->v", n_z_, "0.") << "\n";

    g.comment("Jacobian-times-vector product");
    g << "m->Jv = w; w += " << n_f_ << ";\n";

    g.comment("Unperturbed input");
    g << "m->x0 = w; w += " << n_z_ << ";\n";

    g.comment("Unperturbed output");
    g << "m->f0 = w; w += " << n_f_ << ";\n";

    g.comment("Setup arg, res for calling function");
    g.local("x1", "real_t", "*");
    g << "x1 = m->x;\n";
    for (int j=0; j<n_in; ++j) {
      g << "arg[" << j << "] = x1; x1 += " << derivative_of_.nnz_in(j) << ";\n";
    }
    g.local("f1", "real_t", "*");
    g << "f1 = m->f;\n";
    for (int j=0; j<n_out; ++j) {
      g << "res[" << j << "] = f1; f1 += " << derivative_of_.nnz_out(j) << ";\n";
    }

    g.comment("Loop over derivative directions");
    g.local("i", "int");
    g << "for (i=0; i<" << n_ << "; ++i) {\n";

    g.comment("Copy seeds to v");
    g.local("v1", "real_t", "*");
    g << "v1 = m->v;\n";
    for (int j=0; j<n_in; ++j) {
      int nnz = derivative_of_.nnz_in(j);
      string s = "seed[" + to_string(j) + "]";
      g << "if (" << s << ") " << g.copy(s + "+i*" + to_string(nnz), nnz, "v1") << "\n"
        << "v1 += " << nnz << ";\n";
    }

    g.comment("Reset counter");
    g << "m->n_calls = 0;\n";

    g.comment("Invoke reverse communication algorithm");
    g << "while (" << g.central_diff("m") << ") {\n"
      << "if (" << g(derivative_of_, "arg", "res", "iw", "w") << ") return 1;\n"
      << "}\n";

    g.comment("Gather sensitivities");
    g.local("Jv1", "real_t", "*");
    g << "Jv1 = m->Jv;\n";
    for (int j=0; j<n_out; ++j) {
      int nnz = derivative_of_.nnz_out(j);
      string s = "sens[" + to_string(j) + "]";
      g << "if (" << s << ") " << g.copy("Jv1", nnz, s + "+i*" + to_string(nnz)) << "\n"
        << "Jv1 += " << nnz << ";\n";
    }

    // End for (i=0; ...
    g << "}\n";
  }


} // namespace casadi
