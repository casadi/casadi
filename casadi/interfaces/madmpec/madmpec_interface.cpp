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


#include "madmpec_interface.hpp"
#include <madmpec_runtime_str.h>

#include "casadi/core/casadi_misc.hpp"
#include "../../core/global_options.hpp"
#include "../../core/casadi_interrupt.hpp"
#include "../../core/convexify.hpp"

#include <ctime>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <string>

namespace casadi {

extern "C"
int CASADI_NLPSOL_MADMPEC_EXPORT
casadi_register_nlpsol_madmpec(Nlpsol::Plugin* plugin) {
  plugin->creator = MadmpecInterface::creator;
  plugin->name = "madmpec";
  plugin->doc = MadmpecInterface::meta_doc.c_str();
  plugin->version = CASADI_VERSION;
  plugin->options = &MadmpecInterface::options_;
  plugin->deserialize = &MadmpecInterface::deserialize;
  return 0;
}

extern "C"
void CASADI_NLPSOL_MADMPEC_EXPORT casadi_load_nlpsol_madmpec() {
  Nlpsol::registerPlugin(casadi_register_nlpsol_madmpec);
}

MadmpecInterface::MadmpecInterface(const std::string& name, const Function& nlp)
  : Nlpsol(name, nlp) {
}

MadmpecInterface::~MadmpecInterface() {
  clear_mem();
}

const Options MadmpecInterface::options_
= {{&Nlpsol::options_},
   {{"nw",
     {OT_INTVECTOR,
      "Number of variables"}},
    {"ng",
     {OT_INTVECTOR,
      "Number of constraints"}},
    {"madnlp",
     {OT_DICT,
      "Options to be passed to madnlp"}},
    {"madmpec",
     {OT_DICT,
      "Options to be passed to madmpec"}},
    {"convexify_strategy",
     {OT_STRING,
      "NONE|regularize|eigen-reflect|eigen-clip. "
      "Strategy to convexify the Lagrange Hessian before passing it to the solver."}},
    {"convexify_margin",
     {OT_DOUBLE,
      "When using a convexification strategy, make sure that "
      "the smallest eigenvalue is at least this (default: 1e-7)."}},
    {"ind_cc",
       {OT_INTVECTORVECTOR,
        "List of complementary constraints on simple bounds. "
        "Pair (i, j) encodes complementarity between the bounds on variable i and variable j."}},
    {"cctypes",
       {OT_INTVECTOR,
        "List of complementary constraints on simple bounds. "
        "Pair (i, j) encodes complementarity between the bounds on variable i and variable j."}},
   }
};

void casadi_madmpec_sparsity(const casadi_int* sp, libmad_int *coord_i, libmad_int *coord_j) {
    // convert ccs to cco
    casadi_int ncol = sp[1];
    const casadi_int* colind = sp+2;
    const casadi_int* row = colind+ncol+1;

    for (casadi_int cc=0; cc<ncol; ++cc) {
        for (casadi_int el=colind[cc]; el<colind[cc+1]; ++el) {
            *coord_i++ = row[el]+1;
            *coord_j++ = cc+1;
        }
    }
}

void MadmpecInterface::init(const Dict& opts) {
  // Call the init method of the base class
  Nlpsol::init(opts);
  std::cout << "init" << std::endl;
  casadi_int struct_cnt=0;

  // Default options
  std::string convexify_strategy = "none";
  double convexify_margin = 1e-7;
  casadi_int max_iter_eig = 200;
  std::vector< std::vector<casadi_int> > ind_cc;
  std::vector<casadi_int> cctypes;
  convexify_ = false;

  calc_g_ = true;
  calc_f_ = true;

  // Read options
  for (auto&& op : opts) {
    std::cout << op.first << "  " << op.second << std::endl;
    if (op.first=="convexify_strategy") {
      convexify_strategy = op.second.to_string();
    } else if (op.first=="convexify_margin") {
      convexify_margin = op.second;
    } else if (op.first=="max_iter") {
      max_iter_eig = op.second;
    } else if (op.first=="madmpec") {
      opts_ = op.second;
    } else if (op.first=="madnlpc") {
      mpcc_opts_ = op.second;
    } else if (op.first=="ind_cc") {
      ind_cc = op.second;
    } else if (op.first=="cctypes") {
      cctypes = op.second;
    }
  }
  std::cout << "construct indcc" << std::endl;
  std::cout << ind_cc << std::endl;
  ind_cc1_.reserve(ind_cc.size());
  ind_cc2_.reserve(ind_cc.size());
  cctypes_.reserve(ind_cc.size());
  for (auto && e : ind_cc) {
    casadi_assert(e.size()==2, "Complementary constraints must come in pairs.");
    // TODO(@anton) fix this
    // casadi_assert(e[0]>=1, "Invalid variable index.");
    // casadi_assert(e[1]>=1, "Invalid variable index.");
    // casadi_assert(e[0]<=nx_, "Invalid variable index.");
    // casadi_assert(e[1]<=nx_, "Invalid variable index.");
    ind_cc1_.push_back(e[0]);
    ind_cc2_.push_back(e[1]);
  }
  std::cout << "construct cctypes" << std::endl;
  for (auto && e : cctypes) {
    cctypes_.push_back(e);
  }
  // Do we need second order derivatives?
  exact_hessian_ = true;
  auto hessian_approximation = opts_.find("hessian_approximation");
  if (hessian_approximation!=opts_.end()) {
    exact_hessian_ = hessian_approximation->second == "exact";
  }

  // Setup NLP functions
  std::cout << "construct create_functions" << std::endl;

  create_function("nlp_f", {"x", "p"}, {"f"});
  create_function("nlp_g", {"x", "p"}, {"g"});

  if (!has_function("nlp_grad_f")) {
    create_function("nlp_grad_f", {"x", "p"}, {"grad:f:x"});
  }

  if (!has_function("nlp_jac_g")) {
    create_function("nlp_jac_g", {"x", "p"}, {"jac:g:x"});
  }
  jacg_sp_ = get_function("nlp_jac_g").sparsity_out(0);

  if (!has_function("nlp_hess_l")) {
    create_function("nlp_hess_l", {"x", "p", "lam:f", "lam:g"},
                    {"tril:hess:gamma:x:x"}, {{"gamma", {"f", "g"}}});
  }
  hesslag_sp_ = get_function("nlp_hess_l").sparsity_out(0);
  casadi_assert(hesslag_sp_.is_tril(), "Hessian must be lower triangular");

  if (convexify_strategy!="none") {
    convexify_ = true;
    Dict opts;
    opts["strategy"] = convexify_strategy;
    opts["margin"] = convexify_margin;
    opts["max_iter_eig"] = max_iter_eig;
    opts["verbose"] = verbose_;
    hesslag_sp_ = Convexify::setup(convexify_data_, hesslag_sp_, opts);
  }

  // transform ccs sparsity to cco
  nzj_i_.resize(jacg_sp_.nnz());
  nzj_j_.resize(jacg_sp_.nnz());
  nzh_i_.resize(hesslag_sp_.nnz());
  nzh_j_.resize(hesslag_sp_.nnz());

  casadi_madmpec_sparsity(jacg_sp_, get_ptr(nzj_i_), get_ptr(nzj_j_));
  casadi_madmpec_sparsity(hesslag_sp_, get_ptr(nzh_i_), get_ptr(nzh_j_));

  set_madmpec_prob();

  // Allocate memory
  casadi_int sz_arg, sz_res, sz_w, sz_iw;
  casadi_madmpec_work(&p_, &sz_arg, &sz_res, &sz_iw, &sz_w);

  alloc_arg(sz_arg, true);
  alloc_res(sz_res, true);
  alloc_iw(sz_iw, true);
  alloc_w(sz_w, true);

  std::vector<char*> _argv = {};
  std::string s;

  int argc = _argv.size();
  char** argv = reinterpret_cast<char**>(_argv.data());

}

int MadmpecInterface::init_mem(void* mem) const {
  if (Nlpsol::init_mem(mem)) return 1;
  if (!mem) return 1;
  auto m = static_cast<MadmpecMemory*>(mem);
  std::cout << "init_mem" << std::endl;
  // Now create the new options struct
  std::cout << "creating nlp_opts" << std::endl;
  libmad_create_options_dict(&(m->d.nlp_opts));
  std::cout << "created nlpopts" << std::endl;
  for (const auto& kv : opts_) {
    switch (kv.second.getType()) {
     case OT_DOUBLE:
       libmad_set_double_option(m->d.nlp_opts, kv.first.c_str(), kv.second);
       break;
     case OT_INT:
       libmad_set_int64_option(m->d.nlp_opts, kv.first.c_str(), kv.second.to_int());
       break;
     case OT_STRING:
     {
       std::string s = kv.second.to_string();
       libmad_set_string_option(m->d.nlp_opts, kv.first.c_str(), s.c_str());
     }
     break;
     case OT_BOOL:
       libmad_set_bool_option(m->d.nlp_opts, kv.first.c_str(), kv.second.to_bool());
       break;
     default:
       casadi_error("Unknown option type.");
    }
  }
  std::cout << "creating mpcc_opts" << std::endl;
  // Now create the new options struct
  libmad_create_options_dict(&(m->d.mpcc_opts));
  for (const auto& kv : opts_) {
    switch (kv.second.getType()) {
     case OT_DOUBLE:
       libmad_set_double_option(m->d.mpcc_opts, kv.first.c_str(), kv.second);
       break;
     case OT_INT:
       libmad_set_int64_option(m->d.mpcc_opts, kv.first.c_str(), kv.second.to_int());
       break;
     case OT_STRING:
     {
       std::string s = kv.second.to_string();
       libmad_set_string_option(m->d.mpcc_opts, kv.first.c_str(), s.c_str());
     }
     break;
     case OT_BOOL:
       libmad_set_bool_option(m->d.mpcc_opts, kv.first.c_str(), kv.second.to_bool());
       break;
     default:
       casadi_error("Unknown option type.");
    }
  }
  std::cout << "created mpccopts" << std::endl;

  m->d.ind_cc1 = ind_cc1_.data();
  m->d.ind_cc2 = ind_cc2_.data();
  m->d.cctypes = cctypes_.data();
  m->d.ncc = ind_cc1_.size();
  casadi_madmpec_init_mem(&m->d);

  return 0;
}

void MadmpecInterface::free_mem(void* mem) const {
  auto m = static_cast<MadmpecMemory*>(mem);
  madmpec_free_mem(&m->d);
  delete static_cast<MadmpecMemory*>(mem);
}

/** \brief Set the (persistent) work vectors */
void MadmpecInterface::set_work(void* mem, const double**& arg, double**& res,
                              casadi_int*& iw, double*& w) const {
  auto m = static_cast<MadmpecMemory*>(mem);

  // Set work in base classes
  Nlpsol::set_work(mem, arg, res, iw, w);

  m->d.prob = &p_;
  m->d.nlp = &m->d_nlp;

  casadi_madmpec_init(&m->d, &arg, &res, &iw, &w);

  m->d.nlp->oracle->m = static_cast<void*>(m);
}

int MadmpecInterface::solve(void* mem) const {
  auto m = static_cast<MadmpecMemory*>(mem);

  casadi_madmpec_presolve(&m->d);

  int ret = casadi_madmpec_solve(&m->d);
  if ( ret != 0 ) throw CasadiException("MADMPECError");

  m->success = m->d.success;
  m->unified_return_status = static_cast<UnifiedReturnStatus>(m->d.unified_return_status);

  return 0;
}

Dict MadmpecInterface::get_stats(void* mem) const {
  Dict stats = Nlpsol::get_stats(mem);
  auto m = static_cast<MadmpecMemory*>(mem);
  libmad_int iter, status;
  double primal_feas, dual_feas;
  madnlpc_get_iters(m->d.stats, &iter);
  madnlpc_get_status(m->d.stats, &status);
  madnlpc_get_dual_feas(m->d.stats, &dual_feas);
  madnlpc_get_primal_feas(m->d.stats, &primal_feas);

  stats["iter_count"] = static_cast<casadi_int>(iter);
  Dict madmpec;
  madmpec["dual_feas"] = dual_feas;
  madmpec["primal_feas"] = primal_feas;
  madmpec["status"] = static_cast<casadi_int>(status);
  stats["madmpec"] = madmpec;
  return stats;
}

void MadmpecInterface::set_madmpec_prob() {
  // assign pointer to internal structur casadi_nlp_prob
  // p_nlp_ ~ casadi_nlp_prob casadi internal
  p_.nlp = &p_nlp_;
  // p_ casadi_madmpec_prob

  p_.nnz_jac_g = jacg_sp_.nnz();
  p_.nnz_hess_l = hesslag_sp_.nnz();
  p_.nzj_i = get_ptr(nzj_i_);
  p_.nzj_j = get_ptr(nzj_j_);
  p_.nzh_i = get_ptr(nzh_i_);
  p_.nzh_j = get_ptr(nzh_j_);

  p_.nlp_hess_l = OracleCallback("nlp_hess_l", this);
  p_.nlp_jac_g = OracleCallback("nlp_jac_g", this);
  p_.nlp_grad_f = OracleCallback("nlp_grad_f", this);
  p_.nlp_f = OracleCallback("nlp_f", this);
  p_.nlp_g = OracleCallback("nlp_g", this);

  casadi_madmpec_setup(&p_);
}

void MadmpecInterface::codegen_init_mem(CodeGenerator& g) const {

  g << "libmad_create_options_dict(&(m->d.libmad_opts))";
  for (const auto& kv : opts_) {
    switch (kv.second.getType()) {
     case OT_DOUBLE:
       g << "libmad_set_double_option(" + codegen_mem(g) + ".libmad_opts, " + kv.first + ", " + str(kv.second) + ");\n";
       break;
     case OT_INT:
       g << "libmad_set_int64_option(" + codegen_mem(g) + ".libmad_opts, " + kv.first + ", " + str(kv.second) + ");\n";
       break;
     case OT_STRING:
     {
       std::string s = kv.second.to_string();
       g << "libmad_set_string_option(" + codegen_mem(g) + ".libmad_opts, " + kv.first + ", " + s + ");\n";
     }
     break;
     case OT_BOOL:
       g << "libmad_set_bool_option(" + codegen_mem(g) + ".libmad_opts, " + kv.first + ", " + str(kv.second) + ");\n";
       break;
     default:
       casadi_error("Unknown option type.");
    }
  }
  g << "madmpec_init_mem(&" + codegen_mem(g) + ");\n";
  g << "return 0;\n";
}

void MadmpecInterface::codegen_free_mem(CodeGenerator& g) const {
  // memory deallocation
  g << "madmpec_free_mem(&" + codegen_mem(g) + ");\n";
}

void MadmpecInterface::codegen_declarations(CodeGenerator& g) const {
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
  g.add_include("MadmpecCInterface.h");
}

void MadmpecInterface::codegen_body(CodeGenerator& g) const {
  codegen_body_enter(g);
  g.auxiliaries << g.sanitize_source(madmpec_runtime_str, {"casadi_real"});

  g.local("d", "struct casadi_madmpec_data*");
  g.init_local("d", "&" + codegen_mem(g));
  g.local("p", "struct casadi_madmpec_prob");
  set_madmpec_prob(g);

  g << "casadi_madmpec_init(d, &arg, &res, &iw, &w);\n";
  g << "casadi_oracle_init(d->nlp->oracle, &arg, &res, &iw, &w);\n";
  g << "casadi_madmpec_presolve(d);\n";

  g << "casadi_madmpec_solve(d);\n";

  codegen_body_exit(g);

  if (error_on_fail_) {
    g << "return d->unified_return_status;\n";
  } else {
    g << "return 0;\n";
  }
}

void MadmpecInterface::set_madmpec_prob(CodeGenerator& g) const {
  if (jacg_sp_.size1()>0 && jacg_sp_.nnz()==0) {
    casadi_error("Empty sparsity pattern not supported in MADMPEC C interface");
  }
  g << "d->nlp = &d_nlp;\n";
  g << "d->prob = &p;\n";
  g << "p.nlp = &p_nlp;\n";

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

  g << "casadi_madmpec_setup(&p);\n";
}

MadmpecInterface::MadmpecInterface(DeserializingStream& s) : Nlpsol(s) {
  s.version("MadmpecInterface", 1);
  s.unpack("MadmpecInterface::jacg_sp", jacg_sp_);
  s.unpack("MadmpecInterface::hesslag_sp", hesslag_sp_);
  s.unpack("MadmpecInterface::exact_hessian", exact_hessian_);
  s.unpack("MadmpecInterface::opts", opts_);
  s.unpack("MadmpecInterface::convexify", convexify_);

  s.unpack("MadmpecInterface::nzj_i", nzj_i_);
  s.unpack("MadmpecInterface::nzj_j", nzj_j_);
  s.unpack("MadmpecInterface::nzh_i", nzh_i_);
  s.unpack("MadmpecInterface::nzh_j", nzh_j_);

  set_madmpec_prob();
}

void MadmpecInterface::serialize_body(SerializingStream &s) const {
  Nlpsol::serialize_body(s);
  s.version("MadmpecInterface", 1);

  s.pack("MadmpecInterface::jacg_sp", jacg_sp_);
  s.pack("MadmpecInterface::hesslag_sp", hesslag_sp_);
  s.pack("MadmpecInterface::exact_hessian", exact_hessian_);
  s.pack("MadmpecInterface::opts", opts_);
  s.pack("MadmpecInterface::convexify", convexify_);

  s.pack("MadmpecInterface::nzj_i", nzj_i_);
  s.pack("MadmpecInterface::nzj_j", nzj_j_);
  s.pack("MadmpecInterface::nzh_i", nzh_i_);
  s.pack("MadmpecInterface::nzh_j", nzh_j_);

}

} // namespace casadi
