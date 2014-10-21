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

#include "casadi/core/mx/mx_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/function/mx_function.hpp"

#include "casadi/core/profiling.hpp"
#include "casadi/core/casadi_options.hpp"

#include <iomanip>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_IMPLICITFUNCTION_NEWTON_EXPORT
  casadi_register_implicitfunction_newton(ImplicitFunctionInternal::Plugin* plugin) {
    plugin->creator = Newton::creator;
    plugin->name = "newton";
    plugin->doc = Newton::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_IMPLICITFUNCTION_NEWTON_EXPORT casadi_load_implicitfunction_newton() {
    ImplicitFunctionInternal::registerPlugin(casadi_register_implicitfunction_newton);
  }

  Newton::Newton(const Function& f) : ImplicitFunctionInternal(f) {
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

  void Newton::solveNonLinear() {
    casadi_log("Newton::solveNonLinear:begin");

    // Set up timers for profiling
    double time_zero=0;
    double time_start=0;
    double time_stop=0;
    if (CasadiOptions::profiling && !CasadiOptions::profilingBinary) {
      time_zero = getRealTime();
      CasadiOptions::profilingLog  << "start " << this << ":" <<getOption("name") << std::endl;
    }

    // Pass the inputs to J
    for (int i=0; i<getNumInputs(); ++i) {
      if (i!=iin_) jac_.setInput(input(i), i);
    }

    // Aliases
    DMatrix &u = output(iout_);
    DMatrix &J = jac_.output(0);
    DMatrix &F = jac_.output(1+iout_);

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

      // Print progress
      if (monitored("step") || monitored("stepsize")) {
        std::cout << "Step " << iter << "." << std::endl;
      }

      if (monitored("step")) {
        std::cout << "  u = " << u << std::endl;
      }

      // Use u to evaluate J
      jac_.setInput(u, iin_);
      for (int i=0; i<getNumInputs(); ++i)
        if (i!=iin_) jac_.setInput(input(i), i);

      if (CasadiOptions::profiling) {
        time_start = getRealTime(); // Start timer
      }

      jac_.evaluate();

      // Write out profiling information
      if (CasadiOptions::profiling && !CasadiOptions::profilingBinary) {
        time_stop = getRealTime(); // Stop timer
        CasadiOptions::profilingLog
            << (time_stop-time_start)*1e6 << " ns | "
            << (time_stop-time_zero)*1e3 << " ms | "
            << this << ":" << getOption("name") << ":0|" << jac_.get() << ":"
            << jac_.getOption("name") << "|evaluate jacobian" << std::endl;
      }

      if (monitored("F")) std::cout << "  F = " << F << std::endl;
      if (monitored("normF"))
        std::cout << "  F (min, max, 1-norm, 2-norm) = "
                  << (*std::min_element(F.data().begin(), F.data().end()))
                  << ", " << (*std::max_element(F.data().begin(), F.data().end()))
                  << ", " << sumAll(fabs(F)) << ", " << sqrt(sumAll(F*F)) << std::endl;
      if (monitored("J")) std::cout << "  J = " << J << std::endl;

      double abstol = 0;
      if (numeric_limits<double>::infinity() != abstol_) {
        abstol = std::max((*std::max_element(F.data().begin(),
                                                  F.data().end())),
                               -(*std::min_element(F.data().begin(),
                                                   F.data().end())));
        if (abstol <= abstol_) {
          casadi_log("Converged to acceptable tolerance - abstol: " << abstol_);
          break;
        }
      }

      // Prepare the linear solver with J
      linsol_.setInput(J, LINSOL_A);

      if (CasadiOptions::profiling) {
        time_start = getRealTime(); // Start timer
      }
      linsol_.prepare();
      // Write out profiling information
      if (CasadiOptions::profiling && !CasadiOptions::profilingBinary) {
        time_stop = getRealTime(); // Stop timer
        CasadiOptions::profilingLog
            << (time_stop-time_start)*1e6 << " ns | "
            << (time_stop-time_zero)*1e3 << " ms | "
            << this << ":" << getOption("name")
            << ":1||prepare linear system" << std::endl;
      }

      if (CasadiOptions::profiling) {
        time_start = getRealTime(); // Start timer
      }
      // Solve against F
      linsol_.solve(&F.front(), 1, false);
      if (CasadiOptions::profiling && !CasadiOptions::profilingBinary) {
        time_stop = getRealTime(); // Stop timer
        CasadiOptions::profilingLog
            << (time_stop-time_start)*1e6 << " ns | "
            << (time_stop-time_zero)*1e3 << " ms | "
            << this << ":" << getOption("name") << ":2||solve linear system" << std::endl;
      }

      if (monitored("step")) {
        std::cout << "  step = " << F << std::endl;
      }

      double abstolStep=0;
      if (numeric_limits<double>::infinity() != abstolStep_) {
        abstolStep = std::max((*std::max_element(F.data().begin(),
                                                  F.data().end())),
                               -(*std::min_element(F.data().begin(),
                                                   F.data().end())));
        if (monitored("stepsize")) {
          std::cout << "  stepsize = " << abstolStep << std::endl;
        }
        if (abstolStep <= abstolStep_) {
          casadi_log("Converged to acceptable tolerance - abstolStep: " << abstolStep_);
          break;
        }
      }

      if (print_iteration_) {
        // Only print iteration header once in a while
        if (iter % 10==0) {
          printIteration(std::cout);
        }

        // Print iteration information
        printIteration(std::cout, iter, abstol, abstolStep);
      }

      // Update Xk+1 = Xk - J^(-1) F
      std::transform(u.begin(), u.end(), F.begin(), u.begin(), std::minus<double>());

      // Get auxiliary outputs
      for (int i=0; i<getNumOutputs(); ++i) {
        if (i!=iout_) jac_.getOutput(output(i), 1+i);
      }
    }

    // Store the iteration count
    if (gather_stats_) stats_["iter"] = iter;

    if (success) stats_["return_status"] = "success";

    // Factorization up-to-date
    fact_up_to_date_ = true;

    casadi_log("Newton::solveNonLinear():end after " << iter << " steps");
  }

  void Newton::init() {

    // Call the base class initializer
    ImplicitFunctionInternal::init();

    casadi_assert_message(f_.getNumInputs()>0,
                          "Newton: the supplied f must have at least one input.");
    casadi_assert_message(!linsol_.isNull(),
                          "Newton::init: linear_solver must be supplied");

    if (hasSetOption("max_iter"))
      max_iter_ = getOption("max_iter");

    if (hasSetOption("abstol"))
      abstol_ = getOption("abstol");

    if (hasSetOption("abstolStep"))
      abstolStep_ = getOption("abstolStep");

    print_iteration_ = getOption("print_iteration");

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
