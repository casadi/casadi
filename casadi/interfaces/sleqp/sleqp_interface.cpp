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
#include "sleqp_func.hpp"

// TODO: Pass options
// TODO: Add consistent error handling
// TODO: Take care of exact / quasi Newton
// TODO: Turn Hessian into products?

// TODO: Do we ever need the parameters??


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
    clear_mem();
  }

  const Options SLEQPInterface::options_
  = {{&Nlpsol::options_},
      {{"sleqp",
        {OT_DICT,
        "Options to be passed to SLEQP"}}
      }
  };

  const std::string SLEQPInterface::meta_doc = "";

  void SLEQPInterface::init(const Dict& opts) {
    Nlpsol::init(opts);

    std::cout << "SLEQPInterface::init" << std::endl;

    max_it = SLEQP_NONE;
    time_limit = SLEQP_NONE;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="sleqp") {
        opts_ = op.second;
      }
    }

    auto opt_it = opts.find("sleqp");

    if(opt_it != std::end(opts))
    {
      const Dict& opts_ = opt_it->second;

      for(const auto& opt : opts_)
      {
        if(opt.first == "max_iter")
        {
          max_it = opt.second;

          if(max_it < 0)
          {
            casadi_error("Invalid iteration limit " + std::to_string(max_it));
          }
        }
        else if(opt.first == "max_wall_time")
        {
          time_limit = opt.second;

          if(time_limit < 0.)
          {
            casadi_error("Invalid time limit " + std::to_string(time_limit));
          }
        }
      }

    }

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

    // Allocate work vectors
    alloc_w(nx_, true); // xk_
    alloc_w(ng_, true); // gk_
    alloc_w(nx_, true); // grad_fk_
    alloc_w(jacg_sp_.nnz(), true); // jac_gk_

    if(!fcallback_.is_null())
    {
      // callback xk
      alloc_w(nx_, true);

      // callback lam_xk
      alloc_w(nx_, true);

      // callback lam_gk
      alloc_w(ng_, true);
    }

    // Setup NLP Hessian
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

  void SLEQPInterface::clear_mem_at(SLEQPMemory* m) const {
    std::cout << "SLEQPInterface::clear_mem" << std::endl;

    SLEQP_CALL_EXC(sleqp_solver_release(&m->internal.solver));
    SLEQP_CALL_EXC(sleqp_problem_release(&m->internal.problem));
    SLEQP_CALL_EXC(sleqp_settings_release(&m->internal.settings));
    SLEQP_CALL_EXC(sleqp_vec_free(&m->internal.primal));
  }

  void SLEQPInterface::free_mem(void *mem) const {
    std::cout << "SLEQPInterface::free_mem" << std::endl;

    SLEQPMemory* m = static_cast<SLEQPMemory*>(mem);

    clear_mem_at(m);

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

    check_inputs(m);

    // clear_mem(m);

    casadi_nlpsol_data<double>& d_nlp = m->d_nlp;

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

    // Read options
    for (auto&& op : opts_) {
      bool found = false;
      for (int i = 0; i < SLEQP_NUM_INT_SETTINGS; ++i) {
        SLEQP_SETTINGS_INT ii = static_cast<SLEQP_SETTINGS_INT>(i);
        if (op.first==sleqp_settings_int_name(ii)) {
          found = true;
          SLEQP_CALL_EXC(sleqp_settings_set_int_value(m->internal.settings,
                                                 ii,
                                                 op.second.to_int()));
        }
      }
      for (int i = 0; i < SLEQP_NUM_REAL_SETTINGS; ++i) {
        SLEQP_SETTINGS_REAL ii = static_cast<SLEQP_SETTINGS_REAL>(i);
        if (op.first==sleqp_settings_real_name(ii)) {
          found = true;
          SLEQP_CALL_EXC(sleqp_settings_set_real_value(m->internal.settings,
                                                 ii,
                                                 op.second.to_double()));
        }
      }
      for (int i = 0; i < SLEQP_NUM_BOOL_SETTINGS; ++i) {
        SLEQP_SETTINGS_BOOL ii = static_cast<SLEQP_SETTINGS_BOOL>(i);
        if (op.first==sleqp_settings_bool_name(ii)) {
          found = true;
          SLEQP_CALL_EXC(sleqp_settings_set_bool_value(m->internal.settings,
                                                 ii,
                                                 op.second.to_bool()));
        }
      }
      for (int i = 0; i < SLEQP_NUM_ENUM_SETTINGS; ++i) {
        SLEQP_SETTINGS_ENUM ii = static_cast<SLEQP_SETTINGS_ENUM>(i);
        if (op.first==sleqp_settings_enum_name(ii)) {
          found = true;
          std::string value = op.second.to_string();
          SLEQP_CALL_EXC(sleqp_settings_set_enum_value_by_string(m->internal.settings,
                                                 ii,
                                                 value.c_str()));
        }
      }
      casadi_assert(found, "Could not find option '" + op.first + "'.");
    }

    /*
    SLEQP_CALL_EXC(sleqp_settings_set_enum_value(m->internal.settings,
                                                 SLEQP_SETTINGS_ENUM_DERIV_CHECK,
                                                 SLEQP_DERIV_CHECK_FIRST));
    */

    SLEQP_CALL_EXC(sleqp_problem_create_simple(&m->internal.problem,
                                               func,
                                               var_lb,
                                               var_ub,
                                               cons_lb,
                                               cons_ub,
                                               m->internal.settings));

    SLEQP_CALL_EXC(sleqp_func_release(&func));

    SLEQP_CALL_EXC(sleqp_vec_free(&cons_ub));
    SLEQP_CALL_EXC(sleqp_vec_free(&cons_lb));

    SLEQP_CALL_EXC(sleqp_vec_free(&var_ub));
    SLEQP_CALL_EXC(sleqp_vec_free(&var_lb));


    // No scaling
    SLEQP_CALL_EXC(sleqp_solver_create(&m->internal.solver,
                                       m->internal.problem,
                                       m->internal.primal,
                                       nullptr));

    m->xk = w;
    w += nx_;

    m->gk = w;
    w += ng_;

    m->grad_fk = w;
    w += nx_;

    m->jac_gk = w;
    w += jacg_sp_.nnz();

    if(!fcallback_.is_null())
    {
      m->cb_xk = w;
      w += nx_;

      m->cb_lam_xk = w;
      w += nx_;

      m->cb_lam_gk = w;
      w += ng_;
    }

    /*
    if (exact_hessian_) {
      m->hess_lk = w; w += hesslag_sp_.nnz();
    }
    */

    return;
  }

  UnifiedReturnStatus map_status(SLEQP_STATUS status)
  {
    switch (status) {
    case SLEQP_STATUS_OPTIMAL:
      return SOLVER_RET_SUCCESS;
    case SLEQP_STATUS_INFEASIBLE:
      return SOLVER_RET_INFEASIBLE;

    case SLEQP_STATUS_ABORT_ITER:
    case SLEQP_STATUS_ABORT_MANUAL:
    case SLEQP_STATUS_ABORT_TIME:
      return SOLVER_RET_LIMITED;
    default:
      // case SLEQP_STATUS_UNKNOWN:
      // case SLEQP_STATUS_RUNNING:
      // case SLEQP_STATUS_UNBOUNDED:
      // case SLEQP_STATUS_ABORT_DEADPOINT:
      return SOLVER_RET_UNKNOWN;
    }
  }

  static SLEQP_RETCODE
  accepted_iterate(SleqpSolver* solver, SleqpIterate* iterate, void* data)
  {
    SLEQPMemory* m = static_cast<SLEQPMemory*>(data);

    const SLEQPInterface* interface = m->interface;
    const Function& fcallback_ = interface->fcallback_;

    std::fill_n(m->arg, fcallback_.n_in(), nullptr);

    double ret_double;
    m->res[0] = &ret_double;

    casadi_nlpsol_data<double>& d_nlp = m->d_nlp;

    double obj_val = sleqp_iterate_obj_val(iterate);

    m->arg[NLPSOL_F] = &obj_val;

    SleqpVec* primal = sleqp_iterate_primal(iterate);
    SLEQP_CALL(sleqp_vec_to_raw(primal, m->cb_xk));
    m->arg[NLPSOL_X] = m->cb_xk;

    SleqpVec* cons_val = sleqp_iterate_cons_val(iterate);
    SLEQP_CALL(sleqp_vec_to_raw(cons_val, m->gk));
    m->arg[NLPSOL_G] = m->gk;

    SleqpVec* vars_dual = sleqp_iterate_vars_dual(iterate);
    SLEQP_CALL(sleqp_vec_to_raw(vars_dual, m->cb_lam_xk));
    m->arg[NLPSOL_LAM_X] = m->cb_lam_xk;

    SleqpVec* cons_dual = sleqp_iterate_cons_dual(iterate);
    SLEQP_CALL(sleqp_vec_to_raw(cons_dual, m->cb_lam_gk));
    m->arg[NLPSOL_LAM_G] = m->cb_lam_gk;

    try {
      fcallback_(m->arg, m->res, m->iw, m->w, 0);
    } catch(KeyboardInterruptException& ex) {
      sleqp_raise(SLEQP_CALLBACK_ERROR, "Interrupt caught in callback...");
    } catch(std::exception& ex) {
      casadi_warning("intermediate_callback: " + std::string(ex.what()));
      if (m->iteration_callback_ignore_errors) return SLEQP_OKAY;
      sleqp_raise(SLEQP_CALLBACK_ERROR, "Exception caught in callback...");
    }

    casadi_int ret = static_cast<casadi_int>(ret_double);

    if(ret != 0 && !m->iteration_callback_ignore_errors)
    {
      sleqp_raise(SLEQP_CALLBACK_ERROR, "Error in callback...");
    }

    return SLEQP_OKAY;
  }

  // Solve the NLP
  int SLEQPInterface::solve(void* mem) const {
    std::cout << "SLEQPInterface::solve" << std::endl;

    SLEQPMemory* m = static_cast<SLEQPMemory*>(mem);

    m->iteration_callback_ignore_errors = iteration_callback_ignore_errors_;
    if (!fcallback_.is_null()) {
      SLEQP_CALL_EXC(sleqp_solver_add_callback(m->internal.solver,
                                               SLEQP_SOLVER_EVENT_ACCEPTED_ITERATE,
                                               (void*)accepted_iterate,
                                               mem));
    }

    SLEQP_CALL_EXC(sleqp_solver_solve(m->internal.solver,
                                      max_it,
                                      time_limit));

    SleqpIterate* iterate;

    SLEQP_CALL(sleqp_solver_solution(m->internal.solver, &iterate));

    casadi_nlpsol_data<double>& d_nlp = m->d_nlp;

    m->success = true;
    m->unified_return_status = map_status(sleqp_solver_status(m->internal.solver));

    SleqpVec* primal = sleqp_iterate_primal(iterate);
    SLEQP_CALL_EXC(sleqp_vec_to_raw(primal, d_nlp.z));

    d_nlp.objective = sleqp_iterate_obj_val(iterate);

    SleqpVec* cons_val = sleqp_iterate_cons_val(iterate);
    SLEQP_CALL_EXC(sleqp_vec_to_raw(cons_val, d_nlp.z + nx_));

    SleqpVec* var_dual = sleqp_iterate_vars_dual(iterate);
    SLEQP_CALL_EXC(sleqp_vec_to_raw(var_dual, d_nlp.lam));

    SleqpVec* cons_dual = sleqp_iterate_cons_dual(iterate);
    SLEQP_CALL_EXC(sleqp_vec_to_raw(cons_dual, d_nlp.lam + nx_));

    if (!fcallback_.is_null()) {
      SLEQP_CALL_EXC(sleqp_solver_remove_callback(m->internal.solver,
                                                  SLEQP_SOLVER_EVENT_ACCEPTED_ITERATE,
                                                  (void*)accepted_iterate,
                                                  mem));
    }

    return 0;
  }


  SLEQPInterface::SLEQPInterface(DeserializingStream& s) : Nlpsol(s) {
    s.version("SLEQPInterface", 1);
    s.unpack("SLEQPInterface::jacg_sp", jacg_sp_);
    s.unpack("SLEQPInterface::max_it", max_it);
    s.unpack("SLEQPInterface::time_limit", time_limit);
  }

  void SLEQPInterface::serialize_body(SerializingStream &s) const {
    Nlpsol::serialize_body(s);
    s.version("SLEQPInterface", 1);
    s.pack("SLEQPInterface::jacg_sp", jacg_sp_);
    s.pack("SLEQPInterface::max_it", max_it);
    s.pack("SLEQPInterface::time_limit", time_limit);
  }

} // namespace casadi
