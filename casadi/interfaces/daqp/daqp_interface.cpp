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

#include "daqp_interface.hpp"
#include "casadi/core/nlp_tools.hpp"

#include <daqp_runtime_str.h>
namespace casadi {

  extern "C"
  int CASADI_CONIC_DAQP_EXPORT
  casadi_register_conic_daqp(Conic::Plugin* plugin) {
    plugin->creator = DaqpInterface::creator;
    plugin->name = "daqp";
    plugin->doc = DaqpInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &DaqpInterface::options_;
    plugin->deserialize = &DaqpInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_DAQP_EXPORT casadi_load_conic_daqp() {
    Conic::registerPlugin(casadi_register_conic_daqp);
  }


  DaqpInterface::DaqpInterface(const std::string& name,
                             const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {
  }

  const Options DaqpInterface::options_
  = {{&Conic::options_},
     {{"daqp",
       {OT_DICT,
        "Options to be passed to Daqp."
        }},
     }
   };

  void DaqpInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Conic::init(opts);

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="daqp") {
        opts_ = op.second;
      }
    }

    // Initialize read-only members of class that don't require saving
    // since they can be derived from other read-only members
    set_daqp_prob();

    // Allocate memory
    casadi_int sz_arg, sz_res, sz_w, sz_iw;
    casadi_daqp_work(&p_, &sz_arg, &sz_res, &sz_iw, &sz_w);

    alloc_arg(sz_arg, true);
    alloc_res(sz_res, true);
    alloc_iw(sz_iw, true);
    alloc_w(sz_w, true);
  }

  void DaqpInterface::set_daqp_prob(CodeGenerator& g) const {
    g << "p.qp = &p_qp;\n";
    g << "daqp_default_settings(&p.settings);\n";
    for (auto&& op : opts_) {
      if (op.first=="primal_tol") {
        g << "p.settings.primal_tol = " << op.second.to_double() << ";\n";
      } else if (op.first=="dual_tol") {
        g << "p.settings.dual_tol = " << op.second.to_double() << ";\n";
      } else if (op.first=="zero_tol") {
        g << "p.settings.zero_tol = " << op.second.to_double() << ";\n";
      } else if (op.first=="pivot_tol") {
        g << "p.settings.pivot_tol = " << op.second.to_double() << ";\n";
      } else if (op.first=="progress_tol") {
        g << "p.settings.progress_tol = " << op.second.to_double() << ";\n";
      } else if (op.first=="cycle_tol") {
        g << "p.settings.cycle_tol = " << op.second.to_int() << ";\n";
      } else if (op.first=="iter_limit") {
        g << "p.settings.iter_limit = " << op.second.to_int() << ";\n";
      } else if (op.first=="fval_bound") {
        g << "p.settings.fval_bound = " << op.second.to_double() << ";\n";
      } else if (op.first=="eps_prox") {
        g << "p.settings.eps_prox = " << op.second.to_double() << ";\n";
      } else if (op.first=="eta_prox") {
        g << "p.settings.eta_prox = " << op.second.to_double() << ";\n";
      } else if (op.first=="rho_soft") {
        g << "p.settings.rho_soft = " << op.second.to_double() << ";\n";
      } else if (op.first=="rel_subopt") {
        g << "p.settings.rel_subopt = " << op.second.to_double() << ";\n";
      } else if (op.first=="abs_subopt") {
        g << "p.settings.abs_subopt = " << op.second.to_double() << ";\n";
      } else {
        casadi_error("Unknown option '" + op.first + "'.");
      }
    }

    g << "casadi_daqp_setup(&p);\n";
  }

  void DaqpInterface::codegen_init_mem(CodeGenerator& g) const {
    g << "daqp_init_mem(&" + codegen_mem(g) + ");\n";
    g << "return 0;\n";
  }

  void DaqpInterface::codegen_free_mem(CodeGenerator& g) const {
    g << "daqp_free_mem(&" + codegen_mem(g) + ");\n";
  }

  void DaqpInterface::set_daqp_prob() {
    p_.qp = &p_qp_;

    DAQPSettings* settings = &p_.settings;

    daqp_default_settings(settings);

    for (auto&& op : opts_) {
      if (op.first=="primal_tol") {
        settings->primal_tol = op.second.to_double();
      } else if (op.first=="dual_tol") {
        settings->dual_tol = op.second.to_double();
      } else if (op.first=="zero_tol") {
        settings->zero_tol = op.second.to_double();
      } else if (op.first=="pivot_tol") {
        settings->pivot_tol = op.second.to_double();
      } else if (op.first=="progress_tol") {
        settings->progress_tol = op.second.to_double();
      } else if (op.first=="cycle_tol") {
        settings->cycle_tol = op.second.to_int();
      } else if (op.first=="iter_limit") {
        settings->iter_limit = op.second.to_int();
      } else if (op.first=="fval_bound") {
        settings->fval_bound = op.second.to_double();
      } else if (op.first=="eps_prox") {
        settings->eps_prox = op.second.to_double();
      } else if (op.first=="eta_prox") {
        settings->eta_prox = op.second.to_double();
      } else if (op.first=="rho_soft") {
        settings->rho_soft = op.second.to_double();
      } else if (op.first=="rel_subopt") {
        settings->rel_subopt = op.second.to_double();
      } else if (op.first=="abs_subopt") {
        settings->abs_subopt = op.second.to_double();
      } else {
        casadi_error("Unknown option '" + op.first + "'.");
      }
    }

    casadi_daqp_setup(&p_);
  }

  int DaqpInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    if (!mem) return 1;
    auto m = static_cast<DaqpMemory*>(mem);
    daqp_init_mem(&m->d);

    m->add_stat("preprocessing");
    m->add_stat("solver");
    m->add_stat("postprocessing");

    return 0;
  }

  void DaqpInterface::free_mem(void* mem) const {
    auto m = static_cast<DaqpMemory*>(mem);
    daqp_free_mem(&m->d);
    delete static_cast<DaqpMemory*>(mem);
  }

  /** \brief Set the (persistent) work vectors */
  void DaqpInterface::set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const {

    auto m = static_cast<DaqpMemory*>(mem);

    Conic::set_work(mem, arg, res, iw, w);

    m->d.prob = &p_;
    m->d.qp = &m->d_qp;

    casadi_daqp_init(&m->d, &arg, &res, &iw, &w);

  }

  int DaqpInterface::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<DaqpMemory*>(mem);

    // Statistics
    m->fstats.at("solver").tic();

    casadi_daqp_solve(&m->d, arg, res, iw, w);
    m->fstats.at("solver").toc();

    return 0;
  }

  DaqpInterface::~DaqpInterface() {
    clear_mem();
  }

  void DaqpInterface::codegen_body(CodeGenerator& g) const {
    qp_codegen_body(g);
    g.add_auxiliary(CodeGenerator::AUX_DENSIFY);
    g.add_auxiliary(CodeGenerator::AUX_COPY);
    g.add_include("daqp/api.h");
    g.add_include("stdio.h");

    g.auxiliaries << g.sanitize_source(daqp_runtime_str, {"casadi_real"});

    g.local("d", "struct casadi_daqp_data*");
    g.init_local("d", "&" + codegen_mem(g));
    g.local("p", "struct casadi_daqp_prob");
    set_daqp_prob(g);

    // Setup data structure (corresponds to set_work)
    g << "d->prob = &p;\n";
    g << "d->qp = &d_qp;\n";
    g << "casadi_daqp_init(d, &arg, &res, &iw, &w);\n";

    g << "casadi_daqp_solve(d, arg, res, iw, w);\n";

    g << "if (!d_qp.success) {\n";
    if (error_on_fail_) {
      g << "return -1000;\n";
    } else {
      g << "return -1;\n";
    }
    g << "}\n";
    g << "return 0;\n";
  }

  Dict DaqpInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<DaqpMemory*>(mem);
    stats["return_status"] = m->d.return_status;
    return stats;
  }

  DaqpInterface::DaqpInterface(DeserializingStream& s) : Conic(s) {
    s.version("DaqpInterface", 1);
    s.unpack("DaqpInterface::opts", opts_);
    set_daqp_prob();
  }

  void DaqpInterface::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);

    s.version("DaqpInterface", 1);
    s.pack("DaqpInterface::opts", opts_);
  }

} // end namespace casadi
