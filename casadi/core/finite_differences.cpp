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

  FiniteDiff::FiniteDiff(const std::string& name, int n, double h)
    : FunctionInternal(name), n_(n), h_(h) {
  }

  FiniteDiff::~FiniteDiff() {
  }

  Options FiniteDiff::options_
  = {{&FunctionInternal::options_},
     {{"second_order_stepsize",
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

  void FiniteDiff::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Default options
    h2_ = 1e-3;
    h_max_ = 1.0;
    eps_in_ = eps_out_ = numeric_limits<double>::epsilon();
    u_aim_ = 100;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="second_order_stepsize") {
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
    alloc_res(n_pert(), true); // yk
    alloc_w((n_pert() + 3) * n_y_, true); // yk[:], y0, y, J
    alloc_w(n_z_, true); // z

    // Dimensions
    if (verbose_) {
      casadi_message("Finite differences (" + class_name() + ") with "
                     + str(n_z_) + " inputs, " + str(n_y_)
                     + " outputs and " + str(n_) + " directional derivatives.");
    }

    // Allocate sufficient temporary memory for function evaluation
    alloc(derivative_of_);
  }

  Sparsity FiniteDiff::get_sparsity_in(int i) {
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

  Sparsity FiniteDiff::get_sparsity_out(int i) {
    return repmat(derivative_of_.sparsity_out(i), 1, n_);
  }

  double FiniteDiff::default_in(int ind) const {
    if (ind<derivative_of_.n_in()) {
      return derivative_of_.default_in(ind);
    } else {
      return 0;
    }
  }

  size_t FiniteDiff::get_n_in() {
    return derivative_of_.n_in() + derivative_of_.n_out() + derivative_of_.n_in();
  }

  size_t FiniteDiff::get_n_out() {
    return derivative_of_.n_out();
  }

  std::string FiniteDiff::get_name_in(int i) {
    int n_in = derivative_of_.n_in(), n_out = derivative_of_.n_out();
    if (i<n_in) {
      return derivative_of_.name_in(i);
    } else if (i<n_in+n_out) {
      return "out_" + derivative_of_.name_out(i-n_in);
    } else {
      return "fwd_" + derivative_of_.name_in(i-n_in-n_out);
    }
  }

  std::string FiniteDiff::get_name_out(int i) {
    return "fwd_" + derivative_of_.name_out(i);
  }

  Function CentralDiff::get_forward(int nfwd, const std::string& name,
                                   const std::vector<std::string>& inames,
                                   const std::vector<std::string>& onames,
                                   const Dict& opts) const {
    // Commented out, does not work well
#if 0
    // The second order derivative is calculated as the backwards derivative
    // of the forward derivative, which is equivalent to central differences
    // of second order
    string f_name = "fd_" + name;
    Dict f_opts = {{"derivative_of", derivative_of_}};
    Function f = Function::create(new ForwardDiff(f_name, n_, h_), f_opts);
    // Calculate backwards derivative of f
    f_opts["derivative_of"] = f;
    return Function::create(new ForwardDiff(name, nfwd, -h_), f_opts);
#endif
    return Function::create(new CentralDiff(name, nfwd, h2_), opts);
  }

  int FiniteDiff::eval(const double** arg, double** res, int* iw, double* w, void* mem) const {
    // Shorthands
    int n_in = derivative_of_.n_in(), n_out = derivative_of_.n_out(), n_pert = this->n_pert();

    // Non-differentiated input
    const double** x0 = arg;
    arg += n_in;

    // Non-differentiated output
    double* y0 = w;
    for (int j=0; j<n_out; ++j) {
      const int nnz = derivative_of_.nnz_out(j);
      casadi_copy(*arg++, nnz, w);
      w += nnz;
    }

    // Forward seeds
    const double** seed = arg;
    arg += n_in;

    // Forward sensitivities
    double** sens = res;
    res += n_out;

    // Finite difference approximation
    double* J = w;
    w += n_y_;

    // Perturbed function values
    double** yk = res;
    res += n_pert;
    for (int j=0; j<n_pert; ++j) {
      yk[j] = w, w += n_y_;
    }

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
      for (int k=0; k<n_pert; ++k) {
        // Perturb inputs
        int off = 0;
        for (int j=0; j<n_in; ++j) {
          int nnz = derivative_of_.nnz_in(j);
          casadi_copy(x0[j], nnz, z + off);
          if (seed[j]) casadi_axpy(nnz, pert(k), seed[j] + i*nnz, z + off);
          off += nnz;
        }
        // Evaluate
        if (derivative_of_(arg, res, iw, w)) return 1;
        // Save outputs
        casadi_copy(y, n_y_, yk[k]);
      }

      // Finite difference calculation
      calc_fd(yk, y0, J);

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

  void ForwardDiff::calc_fd(double** yk, double* y0, double* J) const {
    casadi_copy(*yk, n_y_, J);
    casadi_axpy(n_y_, -1., y0, J);
    casadi_scal(n_y_, 1./h_, J);
  }

  void ForwardDiff::calc_fd(CodeGenerator& g, const std::string& yk,
                           const std::string& y0, const std::string& J) const {
    g << g.copy("*"+yk, n_y_, J) << "\n";
    g << g.axpy(n_y_, "-1.", y0, J) << "\n";
    g << g.scal(n_y_, str(1./h_), J) << "\n";
  }

  void CentralDiff::calc_fd(double** yk, double* y0, double* J) const {
    casadi_copy(yk[1], n_y_, J);
    casadi_axpy(n_y_, -1., yk[0], J);
    casadi_scal(n_y_, 0.5/h_, J);
  }

  void CentralDiff::calc_fd(CodeGenerator& g, const std::string& yk,
                           const std::string& y0, const std::string& J) const {
    g << g.copy(yk+"[1]", n_y_, J) << "\n";
    g << g.axpy(n_y_, "-1.", yk+"[0]", J) << "\n";
    g << g.scal(n_y_, str(0.5/h_), J) << "\n";
  }

  void FiniteDiff::codegen_declarations(CodeGenerator& g) const {
    derivative_of_->add_dependency(g);
  }

  void FiniteDiff::codegen_body(CodeGenerator& g) const {
    // Shorthands
    int n_in = derivative_of_.n_in(), n_out = derivative_of_.n_out(), n_pert = this->n_pert();

    g.comment("Non-differentiated input");
    g.local("x0", "const casadi_real", "**");
    g << "x0 = arg, arg += " << n_in << ";\n";

    g.comment("Non-differentiated output");
    g.local("y0", "casadi_real", "*");
    g << "y0 = w;\n";
    for (int j=0; j<n_out; ++j) {
      const int nnz = derivative_of_.nnz_out(j);
      g << g.copy("*arg++", nnz, "w") << " w += " << nnz << ";\n";
    }

    g.comment("Forward seeds");
    g.local("seed", "const casadi_real", "**");
    g << "seed = arg, arg += " << n_in << ";\n";

    g.comment("Forward sensitivities");
    g.local("sens", "casadi_real", "**");
    g << "sens = res, res += " << n_out << ";\n";

    g.comment("Finite difference approximation");
    g.local("J", "casadi_real", "*");
    g << "J = w, w += " << n_y_ << ";\n";

    g.comment("Perturbed function value");
    g.local("yk", "casadi_real", "**");
    g << "yk = res, res += " << n_pert << ";\n";
    g.local("j", "int");
    g << "for (j=0; j<" << n_pert << "; ++j) yk[j] = w, w += " << n_y_ << ";\n";

    g.comment("Setup arg and z for evaluation");
    g.local("z", "casadi_real", "*");
    g << "z = w;\n";
    for (int j=0; j<n_in; ++j) {
      g << "arg[" << j << "] = w, w += " << derivative_of_.nnz_in(j) << ";\n";
    }

    g.comment("Setup res and y for evaluation");
    g.local("y", "casadi_real", "*");
    g << "y = w;\n";
    for (int j=0; j<n_out; ++j) {
      g << "res[" << j << "] = w, w += " << derivative_of_.nnz_out(j) << ";\n";
    }

    g.comment("For all sensitivity directions");
    g.local("i", "int");
    g << "for (i=0; i<" << n_ << "; ++i) {\n";

    g.comment("Calculate perturbed function values");
    g.local("k", "int");
    g << "for (k=0; k<" << n_pert << "; ++k) {\n";

    g.comment("Perturb inputs");
    int off=0;
    for (int j=0; j<n_in; ++j) {
      int nnz = derivative_of_.nnz_in(j);
      string s = "seed[" + str(j) + "]";
      g << g.copy("x0[" + str(j) + "]", nnz, "z+" + str(off)) << "\n"
        << "if ("+s+") " << g.axpy(nnz, pert("k"),
                                   s+"+i*"+str(nnz), "z+" + str(off)) << "\n";
      off += nnz;
    }

    g.comment("Evaluate");
    if (derivative_of_->simplified_call()) {
      g << g(derivative_of_, "z", "y") << ";\n";
    } else {
      g << "if (" << g(derivative_of_, "arg", "res", "iw", "w") << ") return 1;\n";
    }

    g.comment("Save outputs");
    g << g.copy("y", n_y_, "yk[k]") << "\n";

    g << "}\n"; // for (k=0, ...)

    g.comment("Finite difference calculation");
    calc_fd(g, "yk", "y0", "J");

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

  std::string Smoothing::pert(const std::string& k) const {
    string sign = "(2*(" + k + "/2)-1)";
    string len = "(" + k + "%%2+1)";
    return len + "*" + sign + "*" + str(h_);
  }

  double Smoothing::pert(int k) const {
    int sign = 2*(k/2)-1;
    int len = k%2+1;
    return len*sign*h_;
  }

  void Smoothing::calc_fd(double** yk, double* y0, double* J) const {
    // Epsilon to avoid division by zero
    double eps = 1e-14;

    // Split up yk
    double* y_b1 = yk[0]; // one steps back
    double* y_b2 = yk[1]; // two step back
    double* y_f1 = yk[2]; // one steps forward
    double* y_f2 = yk[3]; // two step forward
    // For all components
    for (int i=0; i<n_y_; ++i) {
      // Sum of all non-nan weights
      double sum_weights=0, fd, sm, w;
      J[i] = 0;
      // Forward shifted central differences
      sm = y_f2[i] - 2*y_f1[i] + y0[i];
      if (!isnan(sm)) {
        sum_weights += w = 1./(sm*sm + eps);
        J[i] += w*(y_f2[i]-y0[i]);
      }
      // Central differences
      sm = y_f1[i] - 2*y0[i] + y_b1[i];
      if (!isnan(sm)) {
        sum_weights += w = 2./(sm*sm + eps);
        J[i] += w*(y_f1[i]-y_b1[i]);
      }
      // Backwards shifted central differences
      sm = y0[i] - 2*y_b1[i] + y_b2[i];
      if (!isnan(sm)) {
        sum_weights += w = 1./(sm*sm + eps);
        J[i] += w*(y0[i]-y_b2[i]);
      }
      // Finalize derivative approximation
      J[i] /= sum_weights*2*h_;
    }
  }

  void Smoothing::calc_fd(CodeGenerator& g, const std::string& yk,
                           const std::string& y0, const std::string& J) const {
      casadi_error("not implemented");
//    g << g.copy(yk+"[1]", n_y_, J) << "\n";
  //  g << g.axpy(n_y_, "-1.", yk+"[0]", J) << "\n";
    //g << g.scal(n_y_, str(0.5/h_), J) << "\n";
  }

  Function Smoothing::get_forward(int nfwd, const std::string& name,
                                   const std::vector<std::string>& inames,
                                   const std::vector<std::string>& onames,
                                   const Dict& opts) const {
    return Function::create(new Smoothing(name, nfwd, h2_), opts);
  }

} // namespace casadi
