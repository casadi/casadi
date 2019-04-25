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
    plugin->options = &Newton::options_;
    plugin->deserialize = &Newton::deserialize;
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
    clear_mem();
  }

  const Options Newton::options_
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
        "Print information about each iteration"}},
      {"line_search",
       {OT_BOOL,
        "Enable line-search (default: true)"}}
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
    line_search_ = true;

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
      } else if (op.first=="line_search") {
        line_search_ = op.second;
      }
    }

    casadi_assert(oracle_.n_in()>0,
                          "Newton: the supplied f must have at least one input.");
    casadi_assert(!linsol_.is_null(),
                          "Newton::init: linear_solver must be supplied");

    set_function(oracle_, "g");


    // Allocate memory
    alloc_w(n_, true); // x
    alloc_w(n_, true); // F
    alloc_w(n_, true); // dx trial
    alloc_w(n_, true); // F trial
    alloc_w(sp_jac_.nnz(), true); // J
  }

 void Newton::set_work(void* mem, const double**& arg, double**& res,
                       casadi_int*& iw, double*& w) const {
     Rootfinder::set_work(mem, arg, res, iw, w);
     auto m = static_cast<NewtonMemory*>(mem);
     m->x = w; w += n_;
     m->f = w; w += n_;
     m->x_trial = w; w += n_;
     m->f_trial = w; w += n_;
     m->jac = w; w += sp_jac_.nnz();
  }

  int Newton::solve(void* mem) const {
    auto m = static_cast<NewtonMemory*>(mem);

    scoped_checkout<Linsol> mem_linsol(linsol_);

    // Get the initial guess
    casadi_copy(m->iarg[iin_], n_, m->x);

    // Perform the Newton iterations
    m->iter=0;
    bool success = true;
    while (true) {
      // Break if maximum number of iterations already reached
      if (m->iter >= max_iter_) {
        if (verbose_) casadi_message("Max iterations reached.");
        m->return_status = "max_iteration_reached";
        m->unified_return_status = SOLVER_RET_LIMITED;
        success = false;
        break;
      }

      // Start a new iteration
      m->iter++;

      // Use x to evaluate g and J
      copy_n(m->iarg, n_in_, m->arg);
      m->arg[iin_] = m->x;
      m->res[0] = m->jac;
      copy_n(m->ires, n_out_, m->res+1);
      m->res[1+iout_] = m->f;
      calc_function(m, "jac_f_z");

      // Check convergence
      double abstol = 0;
      if (abstol_ != numeric_limits<double>::infinity()) {
        for (casadi_int i=0; i<n_; ++i) {
          abstol = max(abstol, fabs(m->f[i]));
        }
        if (abstol <= abstol_) {
          if (verbose_) casadi_message("Converged to acceptable tolerance: " + str(abstol_));
          break;
        }
      }

      // Factorize the linear solver with J
      linsol_.nfact(m->jac, mem_linsol);
      linsol_.solve(m->jac, m->f, 1, false, mem_linsol);

      // Check convergence again
      double abstolStep=0;
      if (numeric_limits<double>::infinity() != abstolStep_) {
        for (casadi_int i=0; i<n_; ++i) {
          abstolStep = max(abstolStep, fabs(m->f[i]));
        }
        if (abstolStep <= abstolStep_) {
          if (verbose_) casadi_message("Minimal step size reached: " + str(abstolStep_));
          break;
        }
      }

      double alpha = 1;
      if (line_search_) {
        copy_n(m->iarg, n_in_, m->arg);
        m->arg[iin_] = m->x_trial;
        copy_n(m->ires, n_out_, m->res);
        m->res[iout_] = m->f_trial;
        while (1) {
          // Xtrial = Xk - alpha*J^(-1) F
          copy_n(m->x, n_, m->x_trial);
          casadi_axpy(n_, -alpha, m->f, m->x_trial);
          calc_function(m, "g");

          double abstol_trial = casadi_norm_inf(n_, m->f_trial);
          if (abstol_trial<=(1-alpha/2)*abstol) {
            copy_n(m->x_trial, n_, m->x);
            break;
          }
          if (alpha*abstolStep <= abstolStep_) {
            if (verbose_) casadi_message("Linesearch did not find a descent step "
                                         "for step size " + str(alpha*abstolStep));
            success = false;
            break;
          }
          alpha*= 0.5;
        }
        if (!success) break;
      } else {
        // X = Xk - J^(-1) F
        casadi_axpy(n_, -alpha, m->f, m->x);
      }

      if (print_iteration_) {
        // Only print iteration header once in a while
        if ((m->iter-1) % 10 ==0) {
          printIteration(uout());
        }

        // Print iteration information
        printIteration(uout(), m->iter, abstol, abstolStep, alpha);
      }

    }

    // Get the solution
    casadi_copy(m->x, n_, m->ires[iout_]);

    // Store the iteration count
    if (success) m->return_status = "success";
    if (verbose_) casadi_message("Newton algorithm took " + str(m->iter) + " steps");

    m->success = success;

    return 0;
  }

  void Newton::printIteration(std::ostream &stream) const {
    stream << setw(5) << "iter";
    stream << setw(10) << "res";
    stream << setw(10) << "step";
    if (line_search_) stream << setw(10) << "alpha";
    stream << std::endl;
    stream.unsetf(std::ios::floatfield);
  }

  void Newton::printIteration(std::ostream &stream, casadi_int iter,
                              double abstol, double abstolStep, double alpha) const {

    std::ios_base::fmtflags f = stream.flags();
    stream << setw(5) << iter;
    stream << setw(10) << scientific << setprecision(2) << abstol;
    stream << setw(10) << scientific << setprecision(2) << abstolStep;
    if (line_search_) stream << setw(10) << scientific << setprecision(2) << alpha;

    stream << fixed;
    stream << std::endl;
    stream.flags(f);
  }

  int Newton::init_mem(void* mem) const {
    if (Rootfinder::init_mem(mem)) return 1;
    auto m = static_cast<NewtonMemory*>(mem);
    m->return_status = "";
    m->iter = 0;
    return 0;
  }

  Dict Newton::get_stats(void* mem) const {
    Dict stats = Rootfinder::get_stats(mem);
    auto m = static_cast<NewtonMemory*>(mem);
    stats["return_status"] = m->return_status;
    stats["iter_count"] = m->iter;
    return stats;
  }


  Newton::Newton(DeserializingStream& s) : Rootfinder(s) {
    s.version("Newton", 1);
    s.unpack("Newton::max_iter", max_iter_);
    s.unpack("Newton::abstol", abstol_);
    s.unpack("Newton::abstolStep", abstolStep_);
    s.unpack("Newton::print_iteration", print_iteration_);
    s.unpack("Newton::line_search", line_search_);
  }

  void Newton::serialize_body(SerializingStream &s) const {
    Rootfinder::serialize_body(s);
    s.version("Newton", 1);
    s.pack("Newton::max_iter", max_iter_);
    s.pack("Newton::abstol", abstol_);
    s.pack("Newton::abstolStep", abstolStep_);
    s.pack("Newton::print_iteration", print_iteration_);
    s.pack("Newton::line_search", line_search_);
  }

} // namespace casadi
