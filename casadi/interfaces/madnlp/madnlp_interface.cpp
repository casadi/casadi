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


#include "madnlp_interface.hpp"
#include "casadi/core/casadi_misc.hpp"
#include "../../core/global_options.hpp"
#include "../../core/casadi_interrupt.hpp"
#include "../../core/convexify.hpp"

#include <ctime>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <chrono>

#include "MadnlpCInterface.h"

#include <madnlp_runtime_str.h>

#define CASADI_NLPSOL_MADNLP_EXPORT __attribute__((visibility("default")))


namespace casadi {

extern "C"
int CASADI_NLPSOL_MADNLP_EXPORT
casadi_register_nlpsol_madnlp(Nlpsol::Plugin* plugin) {
  plugin->creator = MadnlpInterface::creator;
  plugin->name = "madnlp";
  plugin->doc = MadnlpInterface::meta_doc.c_str();
  plugin->version = CASADI_VERSION;
  plugin->options = &MadnlpInterface::options_;
  plugin->deserialize = &MadnlpInterface::deserialize;
  return 0;
}

extern "C"
void CASADI_NLPSOL_MADNLP_EXPORT casadi_load_nlpsol_madnlp() {
  Nlpsol::registerPlugin(casadi_register_nlpsol_madnlp);
}

MadnlpInterface::MadnlpInterface(const std::string& name, const Function& nlp)
  : Nlpsol(name, nlp) {
}

MadnlpInterface::~MadnlpInterface() {
  clear_mem();
}

const Options MadnlpInterface::options_ = {
   {&Nlpsol::options_},{
    {"nw",
     {OT_INTVECTOR,
      "Number of variables"}},
    {"ng",
     {OT_INTVECTOR,
      "Number of constraints"}},
    {"madnlp",
     {OT_DICT,
      "Options to be passed to madnlp"}},
    {"convexify_strategy",
     {OT_STRING,
      "NONE|regularize|eigen-reflect|eigen-clip. "
      "Strategy to convexify the Lagrange Hessian before passing it to the solver."}},
    {"convexify_margin",
     {OT_DOUBLE,
      "When using a convexification strategy, make sure that "
      "the smallest eigenvalue is at least this (default: 1e-7)."}},
   }
};

void MadnlpInterface::init(const Dict& opts) {
  // Call the init method of the base class
  Nlpsol::init(opts);

  casadi_int struct_cnt=0;

  // Default options
  std::string convexify_strategy = "none";
  double convexify_margin = 1e-7;
  casadi_int max_iter_eig = 200;
  convexify_ = false;

  calc_g_ = true;
  calc_f_ = true;

  // Read options
  for (auto&& op : opts) {
    if (op.first=="nw") {
      nws_ = op.second;
      struct_cnt++;
    } else if (op.first=="ng") {
      ngs_ = op.second;
      struct_cnt++;
    } else if (op.first=="convexify_strategy") {
      convexify_strategy = op.second.to_string();
    } else if (op.first=="convexify_margin") {
      convexify_margin = op.second;
    } else if (op.first=="max_iter") {
      max_iter_eig = op.second;
    } else if (op.first=="madnlp") {
      opts_ = op.second;
    }
  }

  // Do we need second order derivatives?
  exact_hessian_ = true;
  auto hessian_approximation = opts_.find("hessian_approximation");
  if (hessian_approximation!=opts_.end()) {
    exact_hessian_ = hessian_approximation->second == "exact";
  }

  // Setup NLP functions
  create_function("nlp_f", {"x", "p"}, {"f"});
  create_function("nlp_g", {"x", "p"}, {"g"});

  if (!has_function("nlp_grad_f")) {
    create_function("nlp_grad_f", {"x", "p"}, {"grad:f:x"});
  }
  gradf_sp_ = get_function("nlp_grad_f").sparsity_out(1);

  if (!has_function("nlp_jac_g")) {
    create_function("nlp_jac_g", {"x", "p"}, {"g", "jac:g:x"});
  }
  jacg_sp_ = get_function("nlp_jac_g").sparsity_out(1);

  if (!has_function("nlp_hess_l")) {
    create_function("nlp_hess_l", {"x", "p", "lam:f", "lam:g"},
                    {"grad:gamma:x", "hess:gamma:x:x"}, {{"gamma", {"f", "g"}}});
  }
  hesslag_sp_ = get_function("nlp_hess_l").sparsity_out(1);
  casadi_assert(hesslag_sp_.is_symmetric(), "Hessian must be symmetric");

  if (convexify_strategy!="none") {
    convexify_ = true;
    Dict opts;
    opts["strategy"] = convexify_strategy;
    opts["margin"] = convexify_margin;
    opts["max_iter_eig"] = max_iter_eig;
    opts["verbose"] = verbose_;
    hesslag_sp_ = Convexify::setup(convexify_data_, hesslag_sp_, opts);
  }

  set_madnlp_prob();

  // Allocate memory
  casadi_int sz_arg, sz_res, sz_w, sz_iw;
  casadi_madnlp_work(&p_, &sz_arg, &sz_res, &sz_iw, &sz_w);

  alloc_arg(sz_arg, true);
  alloc_res(sz_res, true);
  alloc_iw(sz_iw, true);
  alloc_w(sz_w, true);
}

int MadnlpInterface::init_mem(void* mem) const {
  if (Nlpsol::init_mem(mem)) return 1;
  if (!mem) return 1;
  auto m = static_cast<MadnlpMemory*>(mem);
  madnlp_init_mem(&m->d);

  return 0;
}

void MadnlpInterface::free_mem(void* mem) const {
  auto m = static_cast<MadnlpMemory*>(mem);
  madnlp_free_mem(&m->d);
  delete static_cast<MadnlpMemory*>(mem);
}

/** \brief Set the (persistent) work vectors */
void MadnlpInterface::set_work(void* mem, const double**& arg, double**& res,
                              casadi_int*& iw, double*& w) const {
  auto m = static_cast<MadnlpMemory*>(mem);

  // Set work in base classes
  Nlpsol::set_work(mem, arg, res, iw, w);

  m->d.prob = &p_;
  m->d.nlp = &m->d_nlp;

  casadi_madnlp_init(&m->d, &arg, &res, &iw, &w);

  m->d.nlp->oracle->m = static_cast<void*>(m);

  // options
}

int MadnlpInterface::solve(void* mem) const {
  auto m = static_cast<MadnlpMemory*>(mem);

  casadi_madnlp_presolve(&m->d);

  for (const auto& kv : opts_) {
    switch (madnlp_c_option_type(kv.first.c_str())) {
      case 0:
        madnlp_c_set_option_double(m->d.solver, kv.first.c_str(), kv.second);
        break;
      case 1:
        madnlp_c_set_option_int(m->d.solver, kv.first.c_str(), kv.second.to_int());
        break;
      case 2:
        madnlp_c_set_option_bool(m->d.solver, kv.first.c_str(), kv.second.to_bool());
        break;
      case 3:
        {
          std::string s = kv.second.to_string();
          madnlp_c_set_option_string(m->d.solver, kv.first.c_str(), s.c_str());
        }
        break;
      case -1:
        casadi_error("Madnlp option not supported: " + kv.first);
      default:
        casadi_error("Unknown option type.");
    }
  }

  casadi_madnlp_solve(&m->d);

  m->success = m->d.success;
  m->unified_return_status = static_cast<UnifiedReturnStatus>(m->d.unified_return_status);

  return 0;
}

Dict MadnlpInterface::get_stats(void* mem) const {
  Dict stats = Nlpsol::get_stats(mem);
  auto m = static_cast<MadnlpMemory*>(mem);
  Dict madnlp;
  madnlp["compute_sd_time"] = m->d.stats.compute_sd_time;
  madnlp["duinf_time"] = m->d.stats.duinf_time;
  madnlp["eval_hess_time"] = m->d.stats.eval_hess_time;
  madnlp["eval_jac_time"] = m->d.stats.eval_jac_time;
  madnlp["eval_cv_time"] = m->d.stats.eval_cv_time;
  madnlp["eval_grad_time"] = m->d.stats.eval_grad_time;
  madnlp["eval_obj_time"] = m->d.stats.eval_obj_time;
  madnlp["initialization_time"] = m->d.stats.initialization_time;
  madnlp["time_total"] = m->d.stats.time_total;
  madnlp["eval_hess_count"] = m->d.stats.eval_hess_count;
  madnlp["eval_jac_count"] = m->d.stats.eval_jac_count;
  madnlp["eval_cv_count"] = m->d.stats.eval_cv_count;
  madnlp["eval_grad_count"] = m->d.stats.eval_grad_count;
  madnlp["eval_obj_count"] = m->d.stats.eval_obj_count;
  madnlp["iterations_count"] = m->d.stats.iterations_count;
  madnlp["return_flag"] = m->d.stats.return_flag;
  stats["madnlp"] = madnlp;
  stats["iter_count"]  =m->d.stats.iterations_count;
  return stats;
}

void MadnlpInterface::set_madnlp_prob() {
  // assign pointer to internal structur casadi_nlp_prob
  // p_nlp_ ~ casadi_nlp_prob casadi internal
  p_.nlp = &p_nlp_;
  // p_ casadi_madnlp_prob

  p_.jac_g_ccs = jacg_sp_;
  p_.hess_l_ccs = hesslag_sp_;
  p_.grad_f_ccs = gradf_sp_;
  get_function("nlp_hess_l").save("nlp_hess_l.casadi");
  p_.nlp_hess_l = OracleCallback("nlp_hess_l", this);
  p_.nlp_jac_g = OracleCallback("nlp_jac_g", this);
  p_.nlp_grad_f = OracleCallback("nlp_grad_f", this);
  p_.nlp_f = OracleCallback("nlp_f", this);
  p_.nlp_g = OracleCallback("nlp_g", this);

  //casadi_madnlp_setup(&p_);
}

void MadnlpInterface::codegen_init_mem(CodeGenerator& g) const {
  g << "madnlp_init_mem(&" + codegen_mem(g) + ");\n";
  g << "return 0;\n";
}

void MadnlpInterface::codegen_free_mem(CodeGenerator& g) const {
  // memory deallocation
  g << "madnlp_free_mem(&" + codegen_mem(g) + ");\n";
}

void MadnlpInterface::codegen_declarations(CodeGenerator& g) const {
  Nlpsol::codegen_declarations(g);
  g.add_auxiliary(CodeGenerator::AUX_NLP);
  g.add_auxiliary(CodeGenerator::AUX_MAX);
  g.add_auxiliary(CodeGenerator::AUX_COPY);
  g.add_auxiliary(CodeGenerator::AUX_PROJECT);
  g.add_auxiliary(CodeGenerator::AUX_SCAL);
  g.add_auxiliary(CodeGenerator::AUX_SPARSITY);
  g.add_auxiliary(CodeGenerator::AUX_ORACLE_CALLBACK);
  g.add_auxiliary(CodeGenerator::AUX_DENSIFY);
  g.add_auxiliary(CodeGenerator::AUX_SPARSIFY);
  g.add_auxiliary(CodeGenerator::AUX_INF);
  g.add_dependency(get_function("nlp_f"));
  g.add_dependency(get_function("nlp_grad_f"));
  g.add_dependency(get_function("nlp_g"));
  g.add_dependency(get_function("nlp_jac_g"));
  g.add_dependency(get_function("nlp_hess_l"));
  g.add_include("MadnlpCInterface.h");
}

void MadnlpInterface::codegen_body(CodeGenerator& g) const {
  codegen_body_enter(g);
  g.auxiliaries << g.sanitize_source(madnlp_runtime_str, {"casadi_real"});

  g.local("d", "struct casadi_madnlp_data*");
  g.init_local("d", "&" + codegen_mem(g));
  g.local("p", "struct casadi_madnlp_prob");
  set_madnlp_prob(g);

  g << "casadi_madnlp_init(d, &arg, &res, &iw, &w);\n";
  g << "casadi_oracle_init(d->nlp->oracle, &arg, &res, &iw, &w);\n";
  g << "casadi_madnlp_presolve(d);\n";

  for (const auto& kv : opts_) {
    switch (madnlp_c_option_type(kv.first.c_str())) {
      case 0:
        g << "madnlp_c_set_option_double(d->solver, \"" + kv.first + "\", "
              + str(kv.second) + ");\n";
        break;
      case 1:
        g << "madnlp_c_set_option_int(d->solver, \"" + kv.first + "\", "
              + str(kv.second.to_int()) + ");\n";
        break;
      case 2:
        g << "madnlp_c_set_option_bool(d->solver, \"" + kv.first + "\", "
              + str(static_cast<int>(kv.second.to_bool())) + ");\n";
        break;
      case 3:
        {
          std::string s = kv.second.to_string();
          g << "madnlp_c_set_option_bool(d->solver, \"" + kv.first + "\", \""
              + s + "\");\n";
        }
        break;
      case -1:
        casadi_error("Madnlp option not supported: " + kv.first);
      default:
        casadi_error("Unknown option type.");
    }
  }

  // Options
  g << "casadi_madnlp_solve(d);\n";

  codegen_body_exit(g);

  if (error_on_fail_) {
    g << "return d->unified_return_status;\n";
  } else {
    g << "return 0;\n";
  }
}

void MadnlpInterface::set_madnlp_prob(CodeGenerator& g) const {
  if (jacg_sp_.size1()>0 && jacg_sp_.nnz()==0) {
    casadi_error("Empty sparsity pattern not supported in MADNLP C interface");
  }
  g << "d->nlp = &d_nlp;\n";
  g << "d->prob = &p;\n";
  g << "p.nlp = &p_nlp;\n";

  g << "p.nw = " << g.constant(nws_) << ";\n";
  g << "p.ng = " << g.constant(ngs_) << ";\n";

  g.setup_callback("p.nlp_jac_g", get_function("nlp_jac_g"));
  g.setup_callback("p.nlp_grad_f", get_function("nlp_grad_f"));
  g.setup_callback("p.nlp_f", get_function("nlp_f"));
  g.setup_callback("p.nlp_g", get_function("nlp_g"));
  g.setup_callback("p.nlp_hess_l", get_function("nlp_hess_l"));

  g << "p.sp_a = " << g.sparsity(jacg_sp_) << ";\n";
  if (exact_hessian_) {
    g << "p.sp_h = " << g.sparsity(hesslag_sp_) << ";\n";
  } else {
    g << "p.sp_h = 0;\n";
  }

  g << "casadi_madnlp_setup(&p);\n";
}

MadnlpInterface::MadnlpInterface(DeserializingStream& s) : Nlpsol(s) {
  s.version("MadnlpInterface", 1);
  s.unpack("MadnlpInterface::jacg_sp", jacg_sp_);
  s.unpack("MadnlpInterface::hesslag_sp", hesslag_sp_);
  s.unpack("MadnlpInterface::exact_hessian", exact_hessian_);
  s.unpack("MadnlpInterface::opts", opts_);
  s.unpack("MadnlpInterface::convexify", convexify_);

  s.unpack("MadnlpInterface::nws", nws_);
  s.unpack("MadnlpInterface::ngs", ngs_);

  set_madnlp_prob();
}

void MadnlpInterface::serialize_body(SerializingStream &s) const {
  Nlpsol::serialize_body(s);
  s.version("MadnlpInterface", 1);

  s.pack("MadnlpInterface::jacg_sp", jacg_sp_);
  s.pack("MadnlpInterface::hesslag_sp", hesslag_sp_);
  s.pack("MadnlpInterface::exact_hessian", exact_hessian_);
  s.pack("MadnlpInterface::opts", opts_);
  s.pack("MadnlpInterface::convexify", convexify_);

  s.pack("MadnlpInterface::nws", nws_);
  s.pack("MadnlpInterface::ngs", ngs_);
}

} // namespace casadi
