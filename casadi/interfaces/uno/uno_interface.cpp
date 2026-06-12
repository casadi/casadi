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

#include "uno_interface.hpp"
#include "casadi/core/casadi_misc.hpp"
#include "casadi/core/casadi_interrupt.hpp"
#include "casadi/core/code_generator.hpp"

#include <uno_runtime_str.h>

#include <cstdlib>
#include <limits>

namespace casadi {

  extern "C"
  int CASADI_NLPSOL_UNO_EXPORT
  casadi_register_nlpsol_uno(Nlpsol::Plugin* plugin) {
    plugin->creator = UnoInterface::creator;
    plugin->name = "uno";
    plugin->doc = UnoInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &UnoInterface::options_;
    plugin->deserialize = &UnoInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_UNO_EXPORT casadi_load_nlpsol_uno() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_uno);
  }

  UnoInterface::UnoInterface(const std::string& name, const Function& nlp)
      : Nlpsol(name, nlp) {}

  UnoInterface::~UnoInterface() {
    clear_mem();
  }

  const Options UnoInterface::options_
  = {{&Nlpsol::options_},
     {{"uno",
       {OT_DICT,
        "Options to be passed to UNO"}}
     }
  };

  UnoMemory::UnoMemory(const UnoInterface& uno_interface)
      : self(uno_interface), NlpsolMemory() {
    return_status = "Unset";
    d_uno.solver = nullptr;
    d_uno.model = nullptr;
  }

  UnoMemory::~UnoMemory() {
    casadi_uno_free_mem<double>(&d_uno);
  }

  void UnoInterface::free_mem(void* mem) const {
    delete static_cast<UnoMemory*>(mem);
  }

  void UnoInterface::init(const Dict& opts) {
    Nlpsol::init(opts);

    calc_f_ = true;
    calc_g_ = true;

    for (auto&& op : opts) {
      if (op.first == "uno") opts_ = op.second;
    }

    create_function("nlp_f", {"x", "p"}, {"f"});
    create_function("nlp_g", {"x", "p"}, {"g"});
    if (!has_function("nlp_grad_f")) {
      create_function("nlp_grad_f", {"x", "p"}, {"grad:f:x"});
    }
    Function gf_jg_fcn = create_function("nlp_jac_g", {"x", "p"}, {"jac:g:x"});
    jacg_sp_ = gf_jg_fcn.sparsity_out(0);

    Function hess_l_fcn = create_function("nlp_hess_l", {"x", "p", "lam:f", "lam:g"},
                                  {"triu:hess:gamma:x:x"},
                                  {{"gamma", {"f", "g"}}});
    hesslag_sp_ = hess_l_fcn.sparsity_out(0);
    casadi_assert(hesslag_sp_.is_triu(), "Hessian must be upper triangular");

    // Hessian-vector product (matrix-free, used by e.g. the LBFGS preset):
    // the forward derivative of grad(gamma) wrt x, seeded along fwd:x.
    Dict final_options;
    final_options["is_diff_in"]  = std::vector<bool>{true, false, false, false};
    final_options["is_diff_out"] = std::vector<bool>{true};
    Dict func_opts;
    func_opts["final_options"] = final_options;
    create_function("nlp_grad_l", {"x", "p", "lam:f", "lam:g"},
                    {"grad:gamma:x"}, {{"gamma", {"f", "g"}}}, func_opts);
    create_forward("nlp_grad_l", 1);  // registers "fwd1_nlp_grad_l"

    {
      auto jr = jacg_sp_.get_row(),    jc = jacg_sp_.get_col();
      auto hr = hesslag_sp_.get_row(), hc = hesslag_sp_.get_col();
      jacobian_row_indices_.assign(jr.begin(), jr.end());
      jacobian_column_indices_.assign(jc.begin(), jc.end());
      hessian_row_indices_.assign(hr.begin(), hr.end());
      hessian_column_indices_.assign(hc.begin(), hc.end());
    }

    placeholder_lb_g_.assign(ng_, -std::numeric_limits<double>::infinity());
    placeholder_ub_g_.assign(ng_,  std::numeric_limits<double>::infinity());

    set_uno_prob();
  }

  void UnoInterface::set_uno_prob() {
    p_uno_.nlp = p_nlp_;
    p_uno_.sp_a = jacg_sp_;
    p_uno_.sp_h = hesslag_sp_;
    p_uno_.jac_row  = jacobian_row_indices_.data();
    p_uno_.jac_col  = jacobian_column_indices_.data();
    p_uno_.n_jac    = static_cast<uno_int>(jacobian_row_indices_.size());
    p_uno_.hess_row = hessian_row_indices_.data();
    p_uno_.hess_col = hessian_column_indices_.data();
    p_uno_.n_hess   = static_cast<uno_int>(hessian_row_indices_.size());
    p_uno_.obj_cb       = &casadi_uno_obj_wrapper<double>;
    p_uno_.obj_grad_cb  = &casadi_uno_obj_grad_wrapper<double>;
    p_uno_.constr_cb    = &casadi_uno_constr_wrapper<double>;
    p_uno_.jac_cb       = &casadi_uno_jac_wrapper<double>;
    p_uno_.hess_cb      = &casadi_uno_hess_wrapper<double>;
    p_uno_.hess_prod_cb = &casadi_uno_hess_prod_wrapper<double>;
    p_uno_.nlp_f         = OracleCallback("nlp_f", this);
    p_uno_.nlp_g         = OracleCallback("nlp_g", this);
    p_uno_.nlp_grad_f    = OracleCallback("nlp_grad_f", this);
    p_uno_.nlp_jac_g     = OracleCallback("nlp_jac_g", this);
    p_uno_.nlp_hess_l    = OracleCallback("nlp_hess_l", this);
    p_uno_.fwd1_nlp_grad_l = OracleCallback("fwd1_nlp_grad_l", this);
  }

  static void set_uno_option(void* solver, const std::string& name, const GenericType& value) {
    if (value.is_bool()) {
      uno_set_solver_bool_option(solver, name.c_str(), value.to_bool());
    } else if (value.is_int()) {
      uno_set_solver_integer_option(solver, name.c_str(), static_cast<uno_int>(value.to_int()));
    } else if (value.is_double()) {
      uno_set_solver_double_option(solver, name.c_str(), value.to_double());
    } else if (value.is_string()) {
      uno_set_solver_string_option(solver, name.c_str(), value.to_string().c_str());
    } else {
      casadi_assert(false, "Unsupported UNO option type for " + name);
    }
  }

  static void insert_casadi_options(void* solver, Dict opts) {
    Dict casadi_options = Options::sanitize(opts);
    // Apply the preset first: it overwrites some options (e.g. it would reset
    // the Hessian model to exact, clobbering an explicit hessian_model=LBFGS).
    for (auto&& op : casadi_options) {
      if (op.first == "preset") {
        uno_set_solver_preset(solver, op.second.to_string().c_str());
        break;
      }
    }
    for (auto&& op : casadi_options) {
      if (op.first != "preset") set_uno_option(solver, op.first, op.second);
    }
  }

  static uno_int casadi_uno_term_cb_cpp(uno_int, uno_int, const double* primals,
      const double* lower_mult, const double* upper_mult, const double* constraint_mult,
      double, double, double, double, void* user_data) {
    constexpr uno_int CONTINUE = 1;
    constexpr uno_int TERMINATE = 0;
    auto* m = static_cast<UnoMemory*>(user_data);  // wired in init_mem
    const UnoInterface& self = m->self;
    if (self.fcallback_.is_null()) return CONTINUE;

    auto* d_nlp = &m->d_nlp;
    casadi_copy(primals, self.nx_, d_nlp->z);
    for (casadi_int i = 0; i < self.nx_; ++i) {
      // Negative sign convention: bound multipliers are added (cf casadi_uno_solve).
      d_nlp->lam[i] = lower_mult[i] + upper_mult[i];
    }
    if (self.ng_ > 0) {
      casadi_copy(constraint_mult, self.ng_, d_nlp->lam + self.nx_);
    }
    std::fill_n(m->arg, self.fcallback_.n_in(), nullptr);
    m->arg[NLPSOL_X]     = d_nlp->z;
    m->arg[NLPSOL_LAM_X] = d_nlp->lam;
    m->arg[NLPSOL_LAM_G] = d_nlp->lam + self.nx_;
    std::fill_n(m->res, self.fcallback_.n_out(), nullptr);
    double ret_double = 0;
    m->res[0] = &ret_double;
    try {
      self.fcallback_(m->arg, m->res, m->iw, m->w, 0);
    } catch (KeyboardInterruptException&) {
      return TERMINATE;
    } catch (std::exception& ex) {
      casadi_warning(std::string("intermediate_callback: ") + ex.what());
      return self.iteration_callback_ignore_errors_ ? CONTINUE : TERMINATE;
    }
    return static_cast<casadi_int>(ret_double) ? TERMINATE : CONTINUE;
  }

  int UnoInterface::init_mem(void* mem) const {
    if (Nlpsol::init_mem(mem)) return 1;
    auto m = static_cast<UnoMemory*>(mem);

    if (verbose_) {
      uno_int uno_major, uno_minor, uno_patch;
      uno_get_version(&uno_major, &uno_minor, &uno_patch);
      casadi_message("Using Uno v" + str(uno_major) + "." + str(uno_minor) + "." + str(uno_patch));
    }

    m->d_uno.prob = &p_uno_;
    casadi_uno_init_mem<double>(&m->d_uno);
    casadi_uno_init_model<double>(&m->d_uno,
        placeholder_lb_g_.data(), placeholder_ub_g_.data());
    // Override the runtime's no-op termination cb with one that fires
    // Nlpsol::fcallback_ each iteration (opti.callback support). Pass UnoMemory*
    // (not casadi_uno_data*) so the cb can reach m->self.fcallback_ etc.
    uno_set_solver_callbacks(m->d_uno.solver, nullptr, &casadi_uno_term_cb_cpp, m);
    insert_casadi_options(m->d_uno.solver, opts_);
    return 0;
  }

  void UnoInterface::set_work(void* mem, const double**& arg, double**& res,
                              casadi_int*& iw, double*& w) const {
    Nlpsol::set_work(mem, arg, res, iw, w);
    auto m = static_cast<UnoMemory*>(mem);
    // Mirror NlpsolMemory's d_nlp into the by-value field on casadi_uno_data
    // (the runtime helpers and codegen path both read d->nlp.* there). Then
    // retarget the prob pointer at our own p_uno_.nlp -- the canonical
    // copy that the codegen path also references.
    m->d_uno.nlp = m->d_nlp;
    m->d_uno.nlp.prob = &p_uno_.nlp;
    // Wiring trap (cf casadi_nlpsol_plugin skill): the OracleCallback path
    // dispatches via cb->oracle_->calc_function(d->m, ...). Without this set,
    // d->m is NULL on the first eval and the plugin segfaults.
    m->d_nlp.oracle->m = static_cast<void*>(m);
    m->d_uno.nlp.oracle = m->d_nlp.oracle;
  }

  inline const char* return_status_string(int sol) {
    switch (sol) {
    case UNO_FEASIBLE_KKT_POINT:          return "Converged with feasible KKT point";
    case UNO_FEASIBLE_FJ_POINT:           return "Converged with feasible FJ point";
    case UNO_INFEASIBLE_STATIONARY_POINT: return "Converged with infeasible stationary point";
    case UNO_FEASIBLE_SMALL_STEP:         return "Terminated with feasible small step";
    case UNO_INFEASIBLE_SMALL_STEP:       return "Terminated with infeasible small step";
    case UNO_UNBOUNDED:                   return "Terminated with unbounded problem";
    case UNO_NOT_OPTIMAL:                 return "Terminated with not optimal point";
    default:                              return "Terminated with an unknown status";
    }
  }

  int UnoInterface::solve(void* mem) const {
    auto m = static_cast<UnoMemory*>(mem);

    casadi_uno_solve<double>(&m->d_uno);

    m->success = m->d_uno.success;
    m->unified_return_status = static_cast<UnifiedReturnStatus>(m->d_uno.unified_return_status);
    m->return_status = return_status_string(uno_get_solution_status(m->d_uno.solver));
    return 0;
  }

  Dict UnoInterface::get_stats(void* mem) const {
    Dict stats = Nlpsol::get_stats(mem);
    auto m = static_cast<UnoMemory*>(mem);
    stats["return_status"]       = m->return_status;
    stats["iter_count"]          = static_cast<casadi_int>(m->d_uno.iter_count);
    stats["primal_infeasbility"] = m->d_uno.primal_infeasibility;
    stats["stationarity"]        = m->d_uno.stationarity;
    stats["complementarity"]     = m->d_uno.complementarity;
    return stats;
  }

  void UnoInterface::serialize_body(SerializingStream& s) const {
    Nlpsol::serialize_body(s);
    s.version("UnoInterface", 1);
    s.pack("UnoInterface::jacg_sp",    jacg_sp_);
    s.pack("UnoInterface::hesslag_sp", hesslag_sp_);
    s.pack("UnoInterface::opts",       opts_);
  }

  UnoInterface::UnoInterface(DeserializingStream& s) : Nlpsol(s) {
    s.version("UnoInterface", 1);
    s.unpack("UnoInterface::jacg_sp",    jacg_sp_);
    s.unpack("UnoInterface::hesslag_sp", hesslag_sp_);
    s.unpack("UnoInterface::opts",       opts_);
    auto jr = jacg_sp_.get_row(),    jc = jacg_sp_.get_col();
    auto hr = hesslag_sp_.get_row(), hc = hesslag_sp_.get_col();
    jacobian_row_indices_.assign(jr.begin(), jr.end());
    jacobian_column_indices_.assign(jc.begin(), jc.end());
    hessian_row_indices_.assign(hr.begin(), hr.end());
    hessian_column_indices_.assign(hc.begin(), hc.end());
    placeholder_lb_g_.assign(ng_, -std::numeric_limits<double>::infinity());
    placeholder_ub_g_.assign(ng_,  std::numeric_limits<double>::infinity());
    set_uno_prob();
  }

  // ----- Codegen --------------------------------------------------------------

  void UnoInterface::codegen_init_mem(CodeGenerator& g) const {
    g.local("d", "struct casadi_uno_data*");
    g.init_local("d", "&" + codegen_mem(g));
    // Static prob: nx/ng + sparsity + callbacks are problem-invariant, so
    // the prob struct lives forever (one per generated function). Storing
    // it function-scope-static lets us call uno_create_model from init_mem,
    // matching the C++ vm path's "build everything at allocation time".
    g.local("p", "static struct casadi_uno_prob");
    set_uno_prob(g);
    g << "d->prob = &p;\n";
    // Wire d->nlp at the persistent NLP scratch on this memory block.
    // (Casadi convention has d_nlp/p_nlp/d_oracle as per-call function-scope
    // locals via Nlpsol::codegen_body_enter; uno opts out of that and uses
    // the by-value fields on casadi_uno_data instead.)
    g << "\n";
    Nlpsol::codegen_setup_constants(g, "d->nlp", "p.nlp", "d->d_oracle");
    g << "casadi_uno_init_mem(d);\n";
    g << "casadi_uno_init_model(d, "
      << g.constant(placeholder_lb_g_) << ", "
      << g.constant(placeholder_ub_g_) << ");\n";
    // Apply user options (statically known at codegen time). The preset must go
    // first: it overwrites some options (e.g. it would reset the Hessian model,
    // clobbering an explicit hessian_model=LBFGS). Mirrors insert_casadi_options.
    for (auto&& kv : opts_) {
      if (kv.first == "preset") {
        g << "uno_set_solver_preset(d->solver, \"" << kv.second.to_string() << "\");\n";
      }
    }
    for (auto&& kv : opts_) {
      const std::string& key = kv.first;
      if (key == "preset") {
        continue;
      } else if (kv.second.is_bool()) {
        g << "uno_set_solver_bool_option(d->solver, \"" << key << "\", "
          << (kv.second.to_bool() ? "1" : "0") << ");\n";
      } else if (kv.second.is_int()) {
        g << "uno_set_solver_integer_option(d->solver, \"" << key << "\", "
          << kv.second.to_int() << ");\n";
      } else if (kv.second.is_double()) {
        g << "uno_set_solver_double_option(d->solver, \"" << key << "\", "
          << g.constant(kv.second.to_double()) << ");\n";
      } else if (kv.second.is_string()) {
        g << "uno_set_solver_string_option(d->solver, \"" << key << "\", \""
          << kv.second.to_string() << "\");\n";
      } else {
        casadi_error("Unsupported uno option type for '" + key + "'");
      }
    }
    g << "return 0;\n";
  }

  void UnoInterface::codegen_free_mem(CodeGenerator& g) const {
    g << "casadi_uno_free_mem(&" + codegen_mem(g) + ");\n";
  }

  void UnoInterface::codegen_declarations(CodeGenerator& g) const {
    Nlpsol::codegen_declarations(g);
    g.add_auxiliary(CodeGenerator::AUX_NLP);
    g.add_auxiliary(CodeGenerator::AUX_ORACLE_CALLBACK);
    g.add_auxiliary(CodeGenerator::AUX_INF);
    g.add_dependency(get_function("nlp_f"));
    g.add_dependency(get_function("nlp_g"));
    g.add_dependency(get_function("nlp_grad_f"));
    g.add_dependency(get_function("nlp_jac_g"));
    g.add_dependency(get_function("nlp_hess_l"));
    g.add_dependency(get_function("fwd1_nlp_grad_l"));
    g.add_include("Uno_C_API.h");
    g.auxiliaries << g.sanitize_source(uno_runtime_str, {"casadi_real"});
  }

  // Point a prob field at a uno_int (32-bit) index array. g.constant only emits
  // casadi_int (64-bit when WITH_LONGLONG_CORE), so emit a typed array via
  // array()+initializer(). It must be "static const": p is a static prob that
  // outlives init_mem, so a plain local array would dangle once init_mem
  // returns. The empty case sets the pointer to 0.
  static void codegen_int_array(CodeGenerator& g, const std::string& name,
      const std::vector<uno_int>& v) {
    if (v.empty()) {
      g << "p." << name << " = 0;\n";
      return;
    }
    g << g.array("static const int", "p_" + name, v.size(),
                 g.initializer(vector_static_cast<casadi_int>(v)));
    g << "p." << name << " = p_" << name << ";\n";
  }

  void UnoInterface::set_uno_prob(CodeGenerator& g) const {
    g << "p.sp_a     = " << g.sparsity(jacg_sp_)    << ";\n";
    g << "p.sp_h     = " << g.sparsity(hesslag_sp_) << ";\n";
    codegen_int_array(g, "jac_row",  jacobian_row_indices_);
    codegen_int_array(g, "jac_col",  jacobian_column_indices_);
    g << "p.n_jac    = " << jacobian_row_indices_.size() << ";\n";
    codegen_int_array(g, "hess_row", hessian_row_indices_);
    codegen_int_array(g, "hess_col", hessian_column_indices_);
    g << "p.n_hess   = " << hessian_row_indices_.size() << ";\n";
    g.setup_callback("p.nlp_f",      get_function("nlp_f"));
    g.setup_callback("p.nlp_g",      get_function("nlp_g"));
    g.setup_callback("p.nlp_grad_f", get_function("nlp_grad_f"));
    g.setup_callback("p.nlp_jac_g",  get_function("nlp_jac_g"));
    g.setup_callback("p.nlp_hess_l", get_function("nlp_hess_l"));
    g.setup_callback("p.fwd1_nlp_grad_l", get_function("fwd1_nlp_grad_l"));
    g << "p.obj_cb       = &casadi_uno_obj_wrapper;\n";
    g << "p.obj_grad_cb  = &casadi_uno_obj_grad_wrapper;\n";
    g << "p.constr_cb    = &casadi_uno_constr_wrapper;\n";
    g << "p.jac_cb       = &casadi_uno_jac_wrapper;\n";
    g << "p.hess_cb      = &casadi_uno_hess_wrapper;\n";
    g << "p.hess_prod_cb = &casadi_uno_hess_prod_wrapper;\n";
  }

  void UnoInterface::codegen_body(CodeGenerator& g) const {
    // No codegen_body_enter / codegen_body_exit: d_nlp / p_nlp / d_oracle
    // live on casadi_uno_data, populated from codegen_init_mem (constants)
    // and the codegen_setup_per_call call below (per-call wiring).
    g.local("d", "struct casadi_uno_data*");
    g.init_local("d", "&" + codegen_mem(g));
    Nlpsol::codegen_setup_per_call(g, "d->nlp");
    g << "casadi_uno_init(d, &arg, &res, &iw, &w);\n";
    g << "casadi_oracle_init(&d->d_oracle, &arg, &res, &iw, &w);\n";
    g << "casadi_uno_solve(d);\n";
    Nlpsol::codegen_post_solve(g, "d->nlp");
    // Signal failure to the caller so error_on_fail works in generated code too
    // (the vm path does this via the return status); otherwise infeasible/limit
    // outcomes look like success.
    if (error_on_fail_) {
      g << "return d->unified_return_status;\n";
    } else {
      g << "return 0;\n";
    }
  }

}  // namespace casadi
