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
        "Differencing scheme [default: 'central']"}},
      {"h_max",
       {OT_DOUBLE,
        "Maximum step size [default 1.0]"}},
      {"eps_in",
       {OT_DOUBLE,
        "Accuracy of function inputs [default: machine precision]"}},
      {"eps_out",
        {OT_DOUBLE,
        "Accuracy of function outputs [default: machine precision]"}},
      {"u_aim",
        {OT_DOUBLE,
        "Target ratio of roundoff error to truncation error [default: 100.]"}}
     }
  };

  void CentralDiff::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Default options
    h_ = 1e-8;
    h2_ = 1e-3;
    h_max_ = 1.0;
    eps_in_ = eps_out_ = numeric_limits<double>::epsilon();
    u_aim_ = 100;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="stepsize") {
        h_ = op.second;
      } else if (op.first=="second_order_stepsize") {
        h2_ = op.second;
      } else if (op.first=="scheme") {
        casadi_warning("Option 'scheme' currently ignored");
      } else if (op.first=="h_max") {
        h_max_ = op.second;
      } else if (op.first=="eps_in") {
        eps_in_ = op.second;
      } else if (op.first=="eps_out") {
        eps_out_ = op.second;
      } else if (op.first=="u_aim") {
        u_aim_ = op.second;
      }
    }

    // Allocate work vector for (perturbed) inputs and outputs
    n_z_ = derivative_of_.nnz_in();
    n_y_ = derivative_of_.nnz_out();
    alloc_w(5 * n_y_, true); // y_pos, y_neg, y0, y, J
    alloc_w(n_z_, true); // z

    // Dimensions
    if (verbose_) {
      casadi_message("Central differences with " + str(n_z_) + " inputs, " + str(n_y_)
                     + " outputs and " + str(n_) + " directional derivatives.");
    }

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

  int CentralDiff::eval(const double** arg, double** res, int* iw, double* w, void* mem) const {
    // Shorthands
    int n_in = derivative_of_.n_in(), n_out = derivative_of_.n_out();

    // Non-differentiated input
    const double** x0 = arg; arg += n_in;

    // Non-differentiated output
    double* y0 = w;
    for (int j=0; j<n_out; ++j) {
      const int nnz = derivative_of_.nnz_out(j);
      casadi_copy(*arg++, nnz, w);
      w += nnz;
    }

    // Forward seeds
    const double** seed = arg; arg += n_in;

    // Forward sensitivities
    double** sens = res; res += n_out;

    // Finite difference approximation
    double* J = w; w += n_y_;

    // Positively perturbed function value
    double* y_pos = w; w += n_y_;

    // Negatively perturned function value
    double* y_neg = w; w += n_y_;

    // Setup arg and z for evaluation
    double *z = w;
    for (int j=0; j<n_in; ++j) {
      arg[j] = w;
      w += derivative_of_.nnz_in(j);
    }

    // Setup res and y for evaluation
    double *y = w;
    for (int j=0; j<n_out; ++j) {
      res[j] = w;
      w += derivative_of_.nnz_out(j);
    }

    // For all sensitivity directions
    for (int i=0; i<n_; ++i) {
      // Calculate perturbed function values
      for (int k=0; k<2; ++k) {
        // Perturb inputs
        int off = 0;
        for (int j=0; j<n_in; ++j) {
          int nnz = derivative_of_.nnz_in(j);
          casadi_copy(x0[j], nnz, z + off);
          if (seed[j]) casadi_axpy(nnz, k ? h_ : -h_, seed[j] + i*nnz, z + off);
          off += nnz;
        }
        // Evaluate
        if (derivative_of_(arg, res, iw, w)) return 1;
        // Save outputs
        casadi_copy(y, n_y_, k ? y_pos : y_neg);
      }

      // Finite difference calculation
      casadi_copy(y_pos, n_y_, J);
      casadi_axpy(n_y_, -1., y_neg, J);
      casadi_scal(n_y_, 0.5/h_, J);

      // Gather sensitivities
      int off = 0;
      for (int j=0; j<n_out; ++j) {
        int nnz = derivative_of_.nnz_out(j);
        if (sens[j]) casadi_copy(J + off, nnz, sens[j] + i*nnz);
        off += nnz;
      }
    }
    return 0;
  }

  void CentralDiff::codegen_declarations(CodeGenerator& g) const {
    derivative_of_->add_dependency(g);
  }

  void CentralDiff::codegen_body(CodeGenerator& g) const {
    // Shorthands
    int n_in = derivative_of_.n_in(), n_out = derivative_of_.n_out();

    g.comment("Non-differentiated input");
    g.local("x0", "const casadi_real", "**");
    g << "x0 = arg; arg += " << n_in << ";\n";

    g.comment("Non-differentiated output");
    g.local("y0", "casadi_real", "*");
    g << "y0 = w;\n";
    for (int j=0; j<n_out; ++j) {
      const int nnz = derivative_of_.nnz_out(j);
      g << g.copy("*arg++", nnz, "w") << " w += " << nnz << ";\n";
    }

    g.comment("Forward seeds");
    g.local("seed", "const casadi_real", "**");
    g << "seed = arg; arg += " << n_in << ";\n";

    g.comment("Forward sensitivities");
    g.local("sens", "casadi_real", "**");
    g << "sens = res; res += " << n_out << ";\n";

    g.comment("Finite difference approximation");
    g.local("J", "casadi_real", "*");
    g << "J = w; w += " << n_y_ << ";\n";

    g.comment("Positively perturbed function value");
    g.local("y_pos", "casadi_real", "*");
    g << "y_pos = w; w += " << n_y_ << ";\n";

    g.comment("Negatively perturned function value");
    g.local("y_neg", "casadi_real", "*");
    g << "y_neg = w; w += " << n_y_ << ";\n";

    g.comment("Setup arg and z for evaluation");
    g.local("z", "casadi_real", "*");
    g << "z = w;\n";
    for (int j=0; j<n_in; ++j) {
      g << "arg[" << j << "] = w; w += " << derivative_of_.nnz_in(j) << ";\n";
    }

    g.comment("Setup res and y for evaluation");
    g.local("y", "casadi_real", "*");
    g << "y = w;\n";
    for (int j=0; j<n_out; ++j) {
      g << "res[" << j << "] = w; w += " << derivative_of_.nnz_out(j) << ";\n";
    }

    g.comment("For all sensitivity directions");
    g.local("i", "int");
    g << "for (i=0; i<" << n_ << "; ++i) {\n";

    g.comment("Calculate perturbed function values");
    g.local("k", "int");
    g << "for (k=0; k<2; ++k) {\n";

    g.comment("Perturb inputs");
    int off=0;
    for (int j=0; j<n_in; ++j) {
      int nnz = derivative_of_.nnz_in(j);
      string s = "seed[" + str(j) + "]";
      string h = "k?" + str(h_) + ":" + str(-h_);
      g << g.copy("x0[" + str(j) + "]", nnz, "z+" + str(off)) << "\n"
        << "if ("+s+") " << g.axpy(nnz, h , s+"+i*"+str(nnz), "z+" + str(off)) << "\n";
      off += nnz;
    }

    g.comment("Evaluate");
    if (derivative_of_->simplified_call()) {
      g << g(derivative_of_, "z", "y") << ";\n";
    } else {
      g << "if (" << g(derivative_of_, "arg", "res", "iw", "w") << ") return 1;\n";
    }

    g.comment("Save outputs");
    g << g.copy("y", n_y_, "k ? y_pos : y_neg") << "\n";

    g << "}\n"; // for (k=0, ...)

    g.comment("Finite difference calculation");
    g << g.copy("y_pos", n_y_, "J") << "\n";
    g << g.axpy(n_y_, "-1.", "y_neg", "J") << "\n";
    g << g.scal(n_y_, str(0.5/h_), "J") << "\n";

    g.comment("Gather sensitivities");
    off = 0;
    for (int j=0; j<n_out; ++j) {
      int nnz = derivative_of_.nnz_out(j);
      string s = "sens[" + str(j) + "]";
      g << "if (" << s << ") " << g.copy("J+" + str(off), nnz, s + "+i*" + str(nnz)) << "\n";
      off += nnz;
    }
    g << "}\n"; // for (i=0, ...)
  }


} // namespace casadi
