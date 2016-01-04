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


#include "newton.hpp"
#include <iomanip>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_ROOTFINDER_NEWTON_EXPORT
  casadi_register_rootfinder_newton(Rootfinder::Plugin* plugin) {
    plugin->creator = Newton::creator;
    plugin->name = "newton";
    plugin->doc = Newton::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_ROOTFINDER_NEWTON_EXPORT casadi_load_rootfinder_newton() {
    Rootfinder::registerPlugin(casadi_register_rootfinder_newton);
  }

  Newton::Newton(const std::string& name, const Function& f)
    : Rootfinder(name, f) {

    addOption("abstol",                      OT_REAL, 1e-12,
              "Stopping criterion tolerance on max(|F|)");
    addOption("abstolStep",                  OT_REAL, 1e-12,
              "Stopping criterion tolerance on step size");
    addOption("max_iter",  OT_INTEGER, 1000,
              "Maximum number of Newton iterations to perform before returning.");
    addOption("monitor",   OT_STRINGVECTOR, GenericType(),  "", "step|stepsize|J|F|normF", true);

    addOption("print_iteration", OT_BOOLEAN, false,
              "Print information about each iteration");
  }

  Newton::~Newton() {
  }

  void Newton::init() {

    // Call the base class initializer
    Rootfinder::init();

    casadi_assert_message(f_.n_in()>0,
                          "Newton: the supplied f must have at least one input.");
    casadi_assert_message(!linsol_.isNull(),
                          "Newton::init: linear_solver must be supplied");

    if (hasSetOption("max_iter"))
      max_iter_ = option("max_iter");

    if (hasSetOption("abstol"))
      abstol_ = option("abstol");

    if (hasSetOption("abstolStep"))
      abstolStep_ = option("abstolStep");

    print_iteration_ = option("print_iteration");

    // Allocate memory
    alloc_w(n_, true); // x
    alloc_w(jac_.nnz_out(1+iout_), true); // F
    alloc_w(jac_.nnz_out(0), true); // J
  }

  Memory* Newton::memory() const {
    return new NewtonMemory();
  }

  void Newton::eval(Memory& mem, const double** arg, double** res,
                    int* iw, double* w) const {
    NewtonMemory& m = dynamic_cast<NewtonMemory&>(mem);

    // IO buffers
    const double** arg1 = arg + n_in();
    double** res1 = res + n_out();

    // Work vectors
    double* x = w; w += n_;
    double* f = w; w += jac_.nnz_out(1+iout_);
    double* jac = w; w += jac_.nnz_out(0);

    // Get the initial guess
    if (arg[iin_]) {
      copy_n(arg[iin_], n_, x);
    } else {
      fill_n(x, n_, 0);
    }

    // Perform the Newton iterations
    int iter=0;
    bool success = true;
    while (true) {
      // Break if maximum number of iterations already reached
      if (iter >= max_iter_) {
        log("eval", "Max. iterations reached.");
        m.return_status = "max_iteration_reached";
        success = false;
        break;
      }

      // Start a new iteration
      iter++;

      // Use x to evaluate J
      copy_n(arg, n_in(), arg1);
      arg1[iin_] = x;
      res1[0] = jac;
      copy_n(res, n_out(), res1+1);
      res1[1+iout_] = f;
      jac_(arg1, res1, iw, w, 0);

      // Check convergence
      double abstol = 0;
      if (abstol_ != numeric_limits<double>::infinity()) {
        for (int i=0; i<n_; ++i) {
          abstol = max(abstol, fabs(f[i]));
        }
        if (abstol <= abstol_) {
          casadi_msg("Converged to acceptable tolerance - abstol: " << abstol_);
          break;
        }
      }

      // Factorize the linear solver with J
      linsol_.setup(arg1 + LINSOL_NUM_IN, res1 + LINSOL_NUM_OUT, iw, w);
      linsol_.linsol_factorize(jac);
      linsol_.linsol_solve(f, 1, false);

      // Check convergence again
      double abstolStep=0;
      if (numeric_limits<double>::infinity() != abstolStep_) {
        for (int i=0; i<n_; ++i) {
          abstolStep = max(abstolStep, fabs(f[i]));
        }
        if (abstolStep <= abstolStep_) {
          casadi_msg("Converged to acceptable tolerance - abstolStep: " << abstolStep_);
          break;
        }
      }

      if (print_iteration_) {
        // Only print iteration header once in a while
        if (iter % 10==0) {
          printIteration(userOut());
        }

        // Print iteration information
        printIteration(userOut(), iter, abstol, abstolStep);
      }

      // Update Xk+1 = Xk - J^(-1) F
      casadi_axpy(n_, -1., f, x);
    }

    // Get the solution
    if (res[iout_]) {
      copy_n(x, n_, res[iout_]);
    }

    // Store the iteration count
    if (gather_stats_) m.iter = iter;
    if (success) m.return_status = "success";

    casadi_msg("Newton::solveNonLinear():end after " << iter << " steps");
  }

  void Newton::printIteration(std::ostream &stream) const {
    stream << setw(5) << "iter";
    stream << setw(10) << "res";
    stream << setw(10) << "step";
    stream << std::endl;
    stream.unsetf(std::ios::floatfield);
  }

  void Newton::printIteration(std::ostream &stream, int iter,
                              double abstol, double abstolStep) const {
    stream << setw(5) << iter;
    stream << setw(10) << scientific << setprecision(2) << abstol;
    stream << setw(10) << scientific << setprecision(2) << abstolStep;

    stream << fixed;
    stream << std::endl;
    stream.unsetf(std::ios::floatfield);
  }

  NewtonMemory::NewtonMemory() {
    return_status = 0;
    iter = 0;
  }

} // namespace casadi
