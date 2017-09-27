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

  FiniteDiff::FiniteDiff(const std::string& name, int n)
    : FunctionInternal(name), n_(n) {
  }

  FiniteDiff::~FiniteDiff() {
  }

  Options FiniteDiff::options_
  = {{&FunctionInternal::options_},
     {{"second_order_stepsize",
       {OT_DOUBLE,
        "Second order perturbation size [default: 1e-3]"}},
      {"h_max",
       {OT_DOUBLE,
        "Maximum step size [default 0]"}},
      {"h_min",
       {OT_DOUBLE,
        "Minimum step size [default inf]"}},
      {"smoothing",
       {OT_DOUBLE,
        "Smoothing regularization [default: machine precision]"}},
      {"reltol",
       {OT_DOUBLE,
        "Accuracy of function inputs [default: query object]"}},
      {"abstol",
        {OT_DOUBLE,
        "Accuracy of function outputs [default: query object]"}},
      {"u_aim",
        {OT_DOUBLE,
        "Target ratio of roundoff error to truncation error [default: 100.]"}},
      {"h_iter",
        {OT_INT,
        "Number of iterations to improve on the step-size "
        "[default: 1 if error estimate available, otherwise 0]"}},
     }
  };

  void FiniteDiff::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Default options
    h_min_ = 0;
    h_max_ = inf;
    smoothing_ = eps;
    reltol_ = derivative_of_->get_reltol();
    abstol_ = derivative_of_->get_abstol();
    h_ = calc_stepsize(abstol_);
    u_aim_ = 100;
    h_iter_ = has_err() ? 1 : 0;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="h") {
        h_ = op.second;
      } else if (op.first=="h_min") {
        h_min_ = op.second;
      } else if (op.first=="h_max") {
        h_max_ = op.second;
      } else if (op.first=="reltol") {
        reltol_ = op.second;
      } else if (op.first=="abstol") {
        abstol_ = op.second;
      } else if (op.first=="smoothing") {
        smoothing_ = op.second;
      } else if (op.first=="u_aim") {
        u_aim_ = op.second;
      } else if (op.first=="h_iter") {
        h_iter_ = op.second;
      }
    }

    // Check h_iter for consistency
    if (h_iter_!=0 && !has_err()) {
      casadi_error("Perturbation size refinement requires an error estimate, "
      "which is not available for the class '" + class_name() + "'. "
      "Choose a different differencing scheme.");
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

  double FiniteDiff::get_default_in(int ind) const {
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
    return Function::create(new CentralDiff(name, nfwd), opts);
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
      // Initial stepsize
      double h = h_;
      // Perform finite difference algorithm with different step sizes
      for (int iter=0; iter<1+h_iter_; ++iter) {
        // Calculate perturbed function values
        for (int k=0; k<n_pert; ++k) {
          // Perturb inputs
          int off = 0;
          for (int j=0; j<n_in; ++j) {
            int nnz = derivative_of_.nnz_in(j);
            casadi_copy(x0[j], nnz, z + off);
            //cout << "k = " << k << ": pert(k, h) = " << pert(k, h) << endl;
            if (seed[j]) casadi_axpy(nnz, pert(k, h), seed[j] + i*nnz, z + off);
            off += nnz;
          }
          // Evaluate
          if (derivative_of_(arg, res, iw, w)) return 1;
          // Save outputs
          casadi_copy(y, n_y_, yk[k]);
        }
        // Finite difference calculation with error estimate
        double u = calc_fd(yk, y0, J, h);
        if (iter==h_iter_) break;

        // Update step size
        if (u < 0) {
          // Perturbation failed, try a smaller step size
          h /= u_aim_;
        } else {
          // Update h to get u near the target ratio
          h *= sqrt(u_aim_ / fmax(1., u));
        }
        // Make sure h stays in the range [h_min_,h_max_]
        h = fmin(fmax(h, h_min_), h_max_);
      }

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

  double ForwardDiff::calc_fd(double** yk, double* y0, double* J, double h) const {
    casadi_copy(*yk, n_y_, J);
    casadi_axpy(n_y_, -1., y0, J);
    casadi_scal(n_y_, 1./h, J);
    return -1;
  }

  void ForwardDiff::calc_fd(CodeGenerator& g, const std::string& yk,
                           const std::string& y0, const std::string& J) const {
    g << g.copy("*"+yk, n_y_, J) << "\n";
    g << g.axpy(n_y_, "-1.", y0, J) << "\n";
    g << g.scal(n_y_, str(1./h_), J) << "\n";
  }

  double CentralDiff::calc_fd(double** yk, double* y0, double* J, double h) const {
    double yf, yc, yb, u=0;
    for (int i=0; i<n_y_; ++i) {
      // Copy to local variables, return -1 if invalid entry
      if (!isfinite((yf=yk[1][i]))) return -1;
      if (!isfinite((yc=y0[i]))) return -1;
      if (!isfinite((yb=yk[0][i]))) return -1;
      // Central difference approximation
      J[i] = (yf - yb)/(2*h);
      // Truncation error
      double err_trunc = yf - 2*yc + yb;
      // Roundoff error
      double err_round = reltol_/h*fmax(fabs(yf-yc), fabs(yc-yb)) + abstol_;
      // Update error estimate
      u = fmax(u, fabs(err_trunc/err_round));
    }
    return u;
  }

  void CentralDiff::calc_fd(CodeGenerator& g, const std::string& yk,
                           const std::string& y0, const std::string& J) const {
    g << g.copy(yk+"[1]", n_y_, J) << "\n";
    g << g.axpy(n_y_, "-1.", yk+"[0]", J) << "\n";
    g << g.scal(n_y_, str(0.5/h_), J) << "\n";
  }

  void FiniteDiff::codegen_declarations(CodeGenerator& g) const {
    g.add_dependency(derivative_of_);
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
    g << "if (" << g(derivative_of_, "arg", "res", "iw", "w") << ") return 1;\n";

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

  double Smoothing::pert(int k, double h) const {
    int sign = 2*(k/2)-1;
    int len = k%2+1;
    return len*sign*h;
  }

  double Smoothing::calc_fd(double** yk, double* y0, double* J, double h) const {
    double err_trunc, err_round, u=0;
    for (int i=0; i<n_y_; ++i) {
      // Reset derivative estimate, sum of weights, error estimate
      double sw=0, ui=0;
      J[i] = 0;
      // Stencil
      double yb, yc, yf;
      // For backward shifted, central and forward shifted
      for (int k=0; k<3; ++k) {
        // Derivative candidate, weight
        double Jk, wk;
        // Calculate shifted finite difference approximation
        if (k==0) {
          // Backward shifted
          // 7.10 in Conte & Carl de Boor: Elementary Numerical Analysis (1972)
          // and 25.3.4 in Abramowitz and Stegun, Handbook of Mathematical Functions (1964)
          if (!isfinite((yc=yk[0][i]))) continue;
          if (!isfinite((yb=yk[1][i]))) continue;
          yf = y0[i];
          Jk = 3*yf - 4*yc + yb;
          wk = 1;
        } else if (k==1) {
          // Central
          // We give this the "nominal weight" 4 since if all weights are equal,
          // this would amount to a five point formula for the derivative
          // (yb2 - 8*yb + 8*yf - y_f2)/(12*h)
          // cf. 25.3.6 in Abramowitz and Stegun, Handbook of Mathematical Functions (1964)
          if (!isfinite((yf=yk[2][i]))) continue;
          if (!isfinite((yb=yc))) continue;
          yc = y0[i];
          Jk = yf-yb;
          wk = 4;
        } else {
          // Forward shifted
          if (!isfinite((yc=yf))) continue;
          if (!isfinite((yf=yk[3][i]))) continue;
          yb = y0[i];
          Jk = -3*yb + 4*yc - yf;
          wk = 1;
        }
        // Truncation error
        double err_trunc = yf - 2*yc + yb;
        // Roundoff error
        double err_round = reltol_/h*fmax(fabs(yf-yc), fabs(yc-yb)) + abstol_;
        // We use the second order derivative as a smoothness measure
        double sm = err_trunc/(h*h);
        // Modify the weight according to smoothness
        wk /= sm*sm + smoothing_;
        sw += wk;
        // Added weighted contribution to weight and error
        J[i] += wk * Jk;
        ui += wk * fabs(err_trunc/err_round);
      }
      // If sw is 0, no stencil worked
      if (sw==0) {
        // Return NaN, ignore weight
        J[i] = nan;
      } else {
        // Finalize estimate using the sum of weights and the step length
        J[i] /= 2*h*sw;
        u = fmax(u, ui/sw);
      }
    }

    return u;
  }

  void Smoothing::calc_fd(CodeGenerator& g, const std::string& yk,
                           const std::string& y0, const std::string& J) const {
      casadi_error("not implemented");
//    g << g.copy(yk+"[1]", n_y_, J) << "\n";
  //  g << g.axpy(n_y_, "-1.", yk+"[0]", J) << "\n";
    //g << g.scal(n_y_, str(0.5/h_), J) << "\n";
  }

  Function ForwardDiff::get_forward(int nfwd, const std::string& name,
                                   const std::vector<std::string>& inames,
                                   const std::vector<std::string>& onames,
                                   const Dict& opts) const {
    return Function::create(new ForwardDiff(name, nfwd), opts);
  }

  Function Smoothing::get_forward(int nfwd, const std::string& name,
                                   const std::vector<std::string>& inames,
                                   const std::vector<std::string>& onames,
                                   const Dict& opts) const {
    return Function::create(new Smoothing(name, nfwd), opts);
  }

} // namespace casadi
