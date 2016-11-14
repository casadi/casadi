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
    plugin->version = CASADI_VERSION;
    return 0;
  }

  extern "C"
  void CASADI_ROOTFINDER_NEWTON_EXPORT casadi_load_rootfinder_newton() {
    Rootfinder::registerPlugin(casadi_register_rootfinder_newton);
  }

  Newton::Newton(const std::string& name, const Function& f)
    : Rootfinder(name, f) {
  }

  Newton::~Newton() {
    clear_memory();
  }

  Options Newton::options_
  = {{&Rootfinder::options_},
     {{"abstol",
       {OT_DOUBLE,
        "Stopping criterion tolerance on max(|F|)"}},
      {"abstolStep",
       {OT_DOUBLE,
        "Stopping criterion tolerance on step size"}},
      {"max_iter",
       {OT_INT,
        "Maximum number of Newton iterations to perform before returning."}},
      {"print_iteration",
       {OT_BOOL,
        "Print information about each iteration"}}
     }
  };

  void Newton::init(const Dict& opts) {

    // Call the base class initializer
    Rootfinder::init(opts);

    // Default options
    max_iter_ = 1000;
    abstol_ = 1e-12;
    abstolStep_ = 1e-12;
    print_iteration_ = false;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="max_iter") {
        max_iter_ = op.second;
      } else if (op.first=="abstol") {
        abstol_ = op.second;
      } else if (op.first=="abstolStep") {
        abstolStep_ = op.second;
      } else if (op.first=="print_iteration") {
        print_iteration_ = op.second;
      }
    }

    casadi_assert_message(oracle_.n_in()>0,
                          "Newton: the supplied f must have at least one input.");
    casadi_assert_message(!linsol_.is_null(),
                          "Newton::init: linear_solver must be supplied");

    // Allocate memory
    alloc_w(n_, true); // x
    alloc_w(n_, true); // F
    alloc_w(sp_jac_.nnz(), true); // J
  }

 void Newton::set_work(void* mem, const double**& arg, double**& res,
                       int*& iw, double*& w) const {
     Rootfinder::set_work(mem, arg, res, iw, w);
     auto m = static_cast<NewtonMemory*>(mem);
     m->x = w; w += n_;
     m->f = w; w += n_;
     m->jac = w; w += sp_jac_.nnz();
  }

  void Newton::solve(void* mem) const {
    auto m = static_cast<NewtonMemory*>(mem);

    // Get the initial guess
    casadi_copy(m->iarg[iin_], n_, m->x);

    // Perform the Newton iterations
    m->iter=0;
    bool success = true;
    while (true) {
      // Break if maximum number of iterations already reached
      if (m->iter >= max_iter_) {
        log("eval", "Max. iterations reached.");
        m->return_status = "max_iteration_reached";
        success = false;
        break;
      }

      // Start a new iteration
      m->iter++;

      // Use x to evaluate J
      copy_n(m->iarg, n_in(), m->arg);
      m->arg[iin_] = m->x;
      m->res[0] = m->jac;
      copy_n(m->ires, n_out(), m->res+1);
      m->res[1+iout_] = m->f;
      calc_function(m, "jac_f_z");

      // Check convergence
      double abstol = 0;
      if (abstol_ != numeric_limits<double>::infinity()) {
        for (int i=0; i<n_; ++i) {
          abstol = max(abstol, fabs(m->f[i]));
        }
        if (abstol <= abstol_) {
          casadi_msg("Converged to acceptable tolerance - abstol: " << abstol_);
          break;
        }
      }

      // Factorize the linear solver with J
      linsol_.factorize(m->jac);
      linsol_.solve(m->f, 1, false);

      // Check convergence again
      double abstolStep=0;
      if (numeric_limits<double>::infinity() != abstolStep_) {
        for (int i=0; i<n_; ++i) {
          abstolStep = max(abstolStep, fabs(m->f[i]));
        }
        if (abstolStep <= abstolStep_) {
          casadi_msg("Converged to acceptable tolerance - abstolStep: " << abstolStep_);
          break;
        }
      }

      if (print_iteration_) {
        // Only print iteration header once in a while
        if (m->iter % 10==0) {
          printIteration(userOut());
        }

        // Print iteration information
        printIteration(userOut(), m->iter, abstol, abstolStep);
      }

      // Update Xk+1 = Xk - J^(-1) F
      casadi_axpy(n_, -1., m->f, m->x);
    }

    // Get the solution
    casadi_copy(m->x, n_, m->ires[iout_]);

    // Store the iteration count
    if (success) m->return_status = "success";

    casadi_msg("Newton::solveNonLinear():end after " << m->iter << " steps");
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

  void Newton::init_memory(void* mem) const {
    Rootfinder::init_memory(mem);
    auto m = static_cast<NewtonMemory*>(mem);
    m->return_status = 0;
    m->iter = 0;
  }

} // namespace casadi
