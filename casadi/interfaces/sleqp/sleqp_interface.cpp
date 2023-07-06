/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#include "sleqp_interface.hpp"
#include "sleqp/pub_problem.h"
#include "sleqp/pub_solver.h"
#include "sleqp/pub_types.h"
#include "sleqp/sparse/pub_vec.h"
#include <cstddef>

// TOOD: Wire up logging functions
// TODO: Convert inf values

// TODO: Use casadi exceptions / error reporting??
#define SLEQP_CALL_EXC(x) \
  do {                                           \
      const SLEQP_RETCODE _status = (x);         \
      if(_status != SLEQP_OKAY) {                \
        throw std::runtime_error("SLEQP error"); \
      }                                          \
  } while(false)


namespace casadi {
  extern "C"
  int CASADI_NLPSOL_SLEQP_EXPORT
  casadi_register_nlpsol_sleqp(Nlpsol::Plugin* plugin) {
    plugin->creator = SLEQPInterface::creator;
    plugin->name = "sleqp";
    plugin->doc = SLEQPInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &SLEQPInterface::options_;
    plugin->deserialize = &SLEQPInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_SLEQP_EXPORT casadi_load_nlpsol_sleqp() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_sleqp);
  }

  SLEQPInterface::SLEQPInterface(const std::string& name, const Function& nlp)
    : Nlpsol(name, nlp) {
  }

  SLEQPInterface::~SLEQPInterface() {
    clear_mem();
  }

  const Options SLEQPInterface::options_
  = {

  };

  const std::string SLEQPInterface::meta_doc = "";

  void SLEQPInterface::init(const Dict& opts) {
    Nlpsol::init(opts);

    std::cout << "SLEQPInterface::init" << std::endl;

    // Setup NLP functions
    create_function("nlp_f", {"x", "p"}, {"f"});
    create_function("nlp_g", {"x", "p"}, {"g"});
    if (!has_function("nlp_grad_f")) {
      create_function("nlp_grad_f", {"x", "p"}, {"f", "grad:f:x"});
    }
    if (!has_function("nlp_jac_g")) {
      create_function("nlp_jac_g", {"x", "p"}, {"g", "jac:g:x"});
    }

    // TODO: Pass options

    // Setup NLP Hessian
    // TODO: Take care of quasi-Newton
    // TODO: Make this a product function
    if (!has_function("nlp_hess_l")) {
      create_function("nlp_hess_l", {"x", "p", "lam:f", "lam:g"},
                      {"triu:hess:gamma:x:x"}, {{"gamma", {"f", "g"}}});
    }
  }

  /** \brief Initalize memory block */
  int SLEQPInterface::init_mem(void* mem) const {
    std::cout << "SLEQPInterface::init_mem" << std::endl;
    if (Nlpsol::init_mem(mem)) return 1;

    SLEQPMemory* m = static_cast<SLEQPMemory*>(mem);

    *m = SLEQPMemory{};

    return 0;
  }

  void SLEQPInterface::clear_mem(SLEQPMemory* m) const {
    SLEQP_CALL_EXC(sleqp_solver_release(&m->solver));
    SLEQP_CALL_EXC(sleqp_problem_release(&m->problem));
    SLEQP_CALL_EXC(sleqp_vec_free(&m->primal));
  }

  void SLEQPInterface::free_mem(void *mem) const {
    SLEQPMemory* m = static_cast<SLEQPMemory*>(mem);

    clear_mem(m);

    delete m;
  }

  /// Get all statistics
  Dict SLEQPInterface::get_stats(void* mem) const {
    std::cout << "SLEQPInterface::get_stats" << std::endl;
    Dict ret = Nlpsol::get_stats(mem);
    return ret;
  }

  /** \brief Set the (persistent) work vectors */
  void SLEQPInterface::set_work(void* mem, const double**& arg, double**& res,
                                casadi_int*& iw, double*& w) const {
    std::cout << "SLEQPInterface::set_work" << std::endl;

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    SLEQPMemory* m = static_cast<SLEQPMemory*>(mem);

    clear_mem(m);

    casadi_nlpsol_data<double> d_nlp = m->d_nlp;

    const int num_vars = nx_;
    const int num_cons = ng_;

    SleqpVec* var_lb;
    SleqpVec* var_ub;

    SLEQP_CALL_EXC(sleqp_vec_create_full(&var_lb, num_vars));
    SLEQP_CALL_EXC(sleqp_vec_create_full(&var_ub, num_vars));

    SLEQP_CALL_EXC(sleqp_vec_set_from_raw(var_lb,
                                          const_cast<double*>(d_nlp.lbx),
                                          num_vars,
                                          0.));

    SLEQP_CALL_EXC(sleqp_vec_set_from_raw(var_ub,
                                          const_cast<double*>(d_nlp.ubx),
                                          num_vars,
                                          0.));

    SLEQP_CALL_EXC(sleqp_vec_create_full(&m->primal, num_vars));

    SLEQP_CALL_EXC(sleqp_vec_set_from_raw(m->primal,
                                          const_cast<double*>(d_nlp.x0),
                                          num_vars,
                                          0.));


    SleqpVec* cons_lb;
    SleqpVec* cons_ub;

    SLEQP_CALL_EXC(sleqp_vec_create_full(&cons_lb, num_cons));
    SLEQP_CALL_EXC(sleqp_vec_create_full(&cons_ub, num_cons));

    SLEQP_CALL_EXC(sleqp_vec_set_from_raw(cons_lb,
                                          const_cast<double*>(d_nlp.lbg),
                                          num_cons,
                                          0.));

    SLEQP_CALL_EXC(sleqp_vec_set_from_raw(cons_ub,
                                          const_cast<double*>(d_nlp.ubg),
                                          num_cons,
                                          0.));

    *m = SLEQPMemory{};

    SleqpFunc* func = NULL;

    SLEQP_CALL_EXC(sleqp_problem_create_simple(&m->problem,
                                               func,
                                               var_lb,
                                               var_ub,
                                               cons_lb,
                                               cons_ub,
                                               m->settings));

    SLEQP_CALL_EXC(sleqp_solver_create(&m->solver, m->problem, m->primal, NULL));

    return;
  }

  // Solve the NLP
  int SLEQPInterface::solve(void* mem) const {
    std::cout << "SLEQPInterface::solve" << std::endl;

    SLEQPMemory* m = static_cast<SLEQPMemory*>(mem);

    // TODO: Pass iteration and time limits
    SLEQP_CALL_EXC(sleqp_solver_solve(m->solver, SLEQP_NONE, SLEQP_NONE));

    //return calc_function(NULL, "nlp_f")==0;

    // TODO: pass result back

    return 0;
  }

} // namespace casadi
