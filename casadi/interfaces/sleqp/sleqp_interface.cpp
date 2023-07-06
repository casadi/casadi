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

#include <cstddef>

#include "sleqp_func.hpp"

// TODO: Pass options
// TODO: Convert inf values

namespace casadi {

  std::string log_level_name(SLEQP_LOG_LEVEL level)
  {
    switch (level) {
      case SLEQP_LOG_DEBUG:
        return "  debug";
      case SLEQP_LOG_INFO:
        return "   info";
      default:
        return "unknown";
    }
  }

  static void casadi_log_output(SLEQP_LOG_LEVEL level,
                                time_t time,
                                const char* message)
  {
    switch(level)
    {
    case SLEQP_LOG_WARN:
      casadi_warning(message);
      break;
    case SLEQP_LOG_ERROR:
      casadi_error(message);
      break;
    default:
      uout() << "[" << log_level_name(level) << "] " << message << std::endl;
    }
  }

  extern "C"
  int CASADI_NLPSOL_SLEQP_EXPORT
  casadi_register_nlpsol_sleqp(Nlpsol::Plugin* plugin) {
    plugin->creator = SLEQPInterface::creator;
    plugin->name = "sleqp";
    plugin->doc = SLEQPInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &SLEQPInterface::options_;
    plugin->deserialize = &SLEQPInterface::deserialize;

    sleqp_log_set_handler(casadi_log_output);
    sleqp_log_set_level(SLEQP_LOG_DEBUG);

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

    jacg_sp_ = get_function("nlp_jac_g").sparsity_out(1);

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

    m->interface = this;

    return 0;
  }

  void SLEQPInterface::clear_mem(SLEQPMemory* m) const {
    std::cout << "SLEQPInterface::clear_mem" << std::endl;

    SLEQP_CALL_EXC(sleqp_solver_release(&m->internal.solver));
    SLEQP_CALL_EXC(sleqp_problem_release(&m->internal.problem));
    SLEQP_CALL_EXC(sleqp_vec_free(&m->internal.primal));

    delete[] m->x;
  }

  void SLEQPInterface::free_mem(void *mem) const {
    std::cout << "SLEQPInterface::free_mem" << std::endl;

    SLEQPMemory* m = static_cast<SLEQPMemory*>(mem);

    clear_mem(m);

    delete m;
  }

  static std::string status_string(SLEQP_STATUS status)
  {
    switch (status) {
    case SLEQP_STATUS_RUNNING:
      return "running";
    case SLEQP_STATUS_OPTIMAL:
      return "optimal";
    case SLEQP_STATUS_INFEASIBLE:
      return "infeasible";
    case SLEQP_STATUS_UNBOUNDED:
      return "unbounded";
    case SLEQP_STATUS_ABORT_DEADPOINT:
      return "deadpoint";
    case SLEQP_STATUS_ABORT_ITER:
      return "iteration limit";
    case SLEQP_STATUS_ABORT_MANUAL:
      return "manual abort";
    case SLEQP_STATUS_ABORT_TIME:
      return "time limit";
    default:
      return "unknown";
    }
  }

  /// Get all statistics
  Dict SLEQPInterface::get_stats(void* mem) const {
    std::cout << "SLEQPInterface::get_stats" << std::endl;

    Dict stats = Nlpsol::get_stats(mem);

    SLEQPMemory* m = static_cast<SLEQPMemory*>(mem);

    SLEQP_STATUS status = sleqp_solver_status(m->internal.solver);

    stats["return_status"] = status_string(status);
    stats["iter_count"] = sleqp_solver_iterations(m->internal.solver);

    return stats;
  }

  /** \brief Set the (persistent) work vectors */
  void SLEQPInterface::set_work(void* mem, const double**& arg, double**& res,
                                casadi_int*& iw, double*& w) const {
    std::cout << "SLEQPInterface::set_work" << std::endl;

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    SLEQPMemory* m = static_cast<SLEQPMemory*>(mem);

    // clear_mem(m);

    casadi_nlpsol_data<double>& d_nlp = m->d_nlp;

    const int num_vars = nx_;
    const int num_cons = ng_;

    m->x = new double[num_vars];

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

    SLEQP_CALL_EXC(sleqp_vec_create_full(&m->internal.primal, num_vars));

    SLEQP_CALL_EXC(sleqp_vec_set_from_raw(m->internal.primal,
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

    SleqpFunc* func = nullptr;

    casadi_sleqp_func_create(&func,
                             num_vars,
                             num_cons,
                             m);

    SLEQP_CALL_EXC(sleqp_settings_create(&m->internal.settings));

    SLEQP_CALL_EXC(sleqp_settings_set_enum_value(m->internal.settings,
                                                 SLEQP_SETTINGS_ENUM_HESS_EVAL,
                                                 SLEQP_HESS_EVAL_DAMPED_BFGS));

    SLEQP_CALL_EXC(sleqp_settings_set_enum_value(m->internal.settings,
                                                 SLEQP_SETTINGS_ENUM_DERIV_CHECK,
                                                 SLEQP_DERIV_CHECK_FIRST));

    SLEQP_CALL_EXC(sleqp_problem_create_simple(&m->internal.problem,
                                               func,
                                               var_lb,
                                               var_ub,
                                               cons_lb,
                                               cons_ub,
                                               m->internal.settings));

    // No scaling
    SLEQP_CALL_EXC(sleqp_solver_create(&m->internal.solver, m->internal.problem, m->internal.primal, nullptr));

    auto jacg_sp_ = get_function("nlp_jac_g").sparsity_out(1);

    m->gk = w;
    w += ng_;

    m->grad_fk = w;
    w += nx_;

    m->jac_gk = w;
    w += jacg_sp_.nnz();

    /*
    if (exact_hessian_) {
      m->hess_lk = w; w += hesslag_sp_.nnz();
    }
    */

    return;
  }

  // Solve the NLP
  int SLEQPInterface::solve(void* mem) const {
    std::cout << "SLEQPInterface::solve" << std::endl;

    SLEQPMemory* m = static_cast<SLEQPMemory*>(mem);

    // TODO: Pass iteration and time limits
    SLEQP_CALL_EXC(sleqp_solver_solve(m->internal.solver, SLEQP_NONE, SLEQP_NONE));

    // TODO: pass result back

    SleqpIterate* iterate;

    SLEQP_CALL(sleqp_solver_solution(m->internal.solver, &iterate));

    casadi_nlpsol_data<double>& d_nlp = m->d_nlp;

    m->success = true;
    m->unified_return_status = SOLVER_RET_SUCCESS;

    SleqpVec* primal = sleqp_iterate_primal(iterate);
    SLEQP_CALL_EXC(sleqp_vec_to_raw(primal, d_nlp.x));

    SLEQP_CALL_EXC(sleqp_vec_to_raw(primal, d_nlp.z));

    d_nlp.objective = sleqp_iterate_obj_val(iterate);
    (*d_nlp.f) = sleqp_iterate_obj_val(iterate);

    SleqpVec* var_dual = sleqp_iterate_vars_dual(iterate);
    SLEQP_CALL_EXC(sleqp_vec_to_raw(var_dual, d_nlp.lam_x));

    SleqpVec* cons_dual = sleqp_iterate_cons_dual(iterate);
    SLEQP_CALL_EXC(sleqp_vec_to_raw(cons_dual, d_nlp.lam_g));

    // TODO: What is lam_p?

    return 0;
  }

} // namespace casadi
