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

#include "casadi/core/profiling.hpp"
#include "casadi/core/casadi_options.hpp"

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

  void Newton::evalD(const double** arg, double** res, int* iw, double* w) {
    casadi_msg("Newton::solveNonLinear:begin");

    // Set up timers for profiling
    double time_zero=0;
    double time_start=0;
    double time_stop=0;
    if (CasadiOptions::profiling && !CasadiOptions::profilingBinary) {
      time_zero = getRealTime();
      CasadiOptions::profilingLog  << "start " << this << ":" <<name_ << std::endl;
    }

    // Get the initial guess
    if (arg[iin_]) {
      copy(arg[iin_], arg[iin_]+nnz_in(iin_), z_.begin());
    } else {
      fill(z_.begin(), z_.end(), 0);
    }

    // Perform the Newton iterations
    int iter=0;

    bool success = true;

    while (true) {
      // Break if maximum number of iterations already reached
      if (iter >= max_iter_) {
        log("evaluate", "Max. iterations reached.");
        stats_["return_status"] = "max_iteration_reached";
        success = false;
        break;
      }

      // Start a new iteration
      iter++;

      // Use u to evaluate J
      jac_.setInputNZ(z_, iin_);

      for (int i=0; i<n_in(); ++i) {
        if (i!=iin_) {
          if (arg[i]) {
            jac_.setInputNZ(arg[i], i);
          } else {
            jac_.setInputNZ(0, i);
          }
        }
      }
      jac_.evaluate();
      DMatrix &J = jac_.output(0);
      DMatrix &F = jac_.output(1+iout_);

      double abstol = 0;
      if (numeric_limits<double>::infinity() != abstol_) {
        abstol = std::max((*std::max_element(F.data().begin(),
                                             F.data().end())),
                          -(*std::min_element(F.data().begin(),
                                              F.data().end())));
        if (abstol <= abstol_) {
          casadi_msg("Converged to acceptable tolerance - abstol: " << abstol_);
          break;
        }
      }

      // Prepare the linear solver with J
      linsol_.setInput(J, LINSOL_A);
      linsol_.prepare();
      linsol_.solve(&F.front(), 1, false);

      double abstolStep=0;
      if (numeric_limits<double>::infinity() != abstolStep_) {
        abstolStep = std::max((*std::max_element(F.data().begin(),
                                                 F.data().end())),
                              -(*std::min_element(F.data().begin(),
                                                  F.data().end())));
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
      std::transform(z_.begin(), z_.end(), F.begin(), z_.begin(), std::minus<double>());
    }

    // Get the solution
    if (res[iout_]) {
      copy_n(z_.begin(), nnz_out(iout_), res[iout_]);
    }

    // Get auxiliary outputs
    for (int i=0; i<n_out(); ++i) {
      if (i!=iout_ && res[i]) {
        copy_n(jac_.output(i+1).ptr(), nnz_out(i), res[i]);
      }
    }

    // Store the iteration count
    if (gather_stats_) stats_["iter"] = iter;
    if (success) stats_["return_status"] = "success";

    casadi_msg("Newton::solveNonLinear():end after " << iter << " steps");
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
    z_.resize(n_);
  }

  void Newton::printIteration(std::ostream &stream) {
    stream << setw(5) << "iter";
    stream << setw(10) << "res";
    stream << setw(10) << "step";
    stream << std::endl;
    stream.unsetf(std::ios::floatfield);
  }

  void Newton::printIteration(std::ostream &stream, int iter, double abstol, double abstolStep) {
    stream << setw(5) << iter;
    stream << setw(10) << scientific << setprecision(2) << abstol;
    stream << setw(10) << scientific << setprecision(2) << abstolStep;

    stream << fixed;
    stream << std::endl;
    stream.unsetf(std::ios::floatfield);
  }


} // namespace casadi
