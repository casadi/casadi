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

#include "highs_interface.hpp"
#include "casadi/core/nlp_tools.hpp"

#include "highs_runtime_str.h"
namespace casadi {

  using namespace std;

  extern "C"
  int CASADI_CONIC_HIGHS_EXPORT
  casadi_register_conic_highs(Conic::Plugin* plugin) {
    plugin->creator = HighsInterface::creator;
    plugin->name = "highs";
    plugin->doc = HighsInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &HighsInterface::options_;
    plugin->deserialize = &HighsInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_HIGHS_EXPORT casadi_load_conic_highs() {
    Conic::registerPlugin(casadi_register_conic_highs);
  }


  HighsInterface::HighsInterface(const std::string& name,
                             const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {
  }

  const Options HighsInterface::options_
  = {{&Conic::options_},
     {{"highs",
       {OT_DICT,
        "Options to be passed to HiGHS."
        }},
     }
   };

  // Options with description https://www.maths.ed.ac.uk/hall/HiGHS/HighsOptions.html
  std::list<std::string> HighsInterface::param_bool = {
    "output_flag",
    "log_to_console",
    "write_solution_to_file",
    "mip_detect_symmetry"
  };

  void HighsInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Conic::init(opts);

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="highs") {
        opts_ = op.second;
      }
    }

    // Initialize read-only members of class that don't require saving
    // since they can be derived from other read-only members
    init_dependent();
    set_highs_prob();

    // Allocate memory
    casadi_int sz_arg, sz_res, sz_w, sz_iw;
    casadi_highs_work(&p_, &sz_arg, &sz_res, &sz_iw, &sz_w);

    alloc_arg(sz_arg, true);
    alloc_res(sz_res, true);
    alloc_iw(sz_iw, true);
    alloc_w(sz_w, true);
  }

  void HighsInterface::init_dependent() {
    colinda_.resize(A_.size2()+1);
    rowa_.resize(A_.nnz());
    colindh_.resize(H_.size2()+1);
    rowh_.resize(H_.nnz());
    copy_vector(A_.colind(), colinda_);
    copy_vector(A_.row(), rowa_);
    copy_vector(H_.colind(), colindh_);
    copy_vector(H_.row(), rowh_);
    if (!discrete_.empty()) {
      integrality_.resize(nx_);
      assign_vector(discrete_, integrality_);
    }
  }


  void codegen_local(CodeGenerator& g, const std::string& name, const std::vector<int>& v) {
    std::string n = name + "[]";
    g.local(n, "static const int");
    std::stringstream init;
    init << "{";
    for (casadi_int i=0;i<v.size();++i) {
      init << v[i];
      if (i<v.size()-1) init << ", ";
    }
    // ISO C forbids empty initializer braces
    if (v.empty()) init << "0";
    init << "}";
    g.init_local(n, init.str());
  }

  void HighsInterface::set_highs_prob(CodeGenerator& g) const {
    g << "p.qp = &p_qp;\n";
    codegen_local(g, "colinda", colinda_);
    codegen_local(g, "rowa", rowa_);
    codegen_local(g, "colindh", colindh_);
    codegen_local(g, "rowh", rowh_);
    if (!discrete_.empty()) {
      codegen_local(g, "integrality", integrality_);
    }
    g << "p.colinda = colinda;\n";
    g << "p.rowa = rowa;\n";
    g << "p.colindh = colindh;\n";
    g << "p.rowh = rowh;\n";
    if (discrete_.empty()) {
      g << "p.integrality = 0;\n";
    } else {
      g << "p.integrality = integrality;\n";
    }
    g << "casadi_highs_setup(&p);\n";
  }

  void HighsInterface::codegen_init_mem(CodeGenerator& g) const {
    g << "highs_init_mem(&" + codegen_mem(g) + ");\n";
    g << "return 0;\n";
  }

  void HighsInterface::codegen_free_mem(CodeGenerator& g) const {
    g << "highs_free_mem(&" + codegen_mem(g) + ");\n";
  }

  void HighsInterface::set_highs_prob() {
    p_.qp = &p_qp_;
    p_.colinda  = get_ptr(colinda_);
    p_.rowa  = get_ptr(rowa_);
    p_.colindh  = get_ptr(colindh_);
    p_.rowh  = get_ptr(rowh_);
    p_.integrality  = get_ptr(integrality_);

    casadi_highs_setup(&p_);
  }

  int HighsInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    if (!mem) return 1;
    auto m = static_cast<HighsMemory*>(mem);
    highs_init_mem(&m->d);

    m->add_stat("preprocessing");
    m->add_stat("solver");
    m->add_stat("postprocessing");

    return 0;
  }

  void HighsInterface::free_mem(void* mem) const {
    auto m = static_cast<HighsMemory*>(mem);
    highs_free_mem(&m->d);
    delete static_cast<HighsMemory*>(mem);
   }

  /** \brief Set the (persistent) work vectors */
  void HighsInterface::set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const {

    auto m = static_cast<HighsMemory*>(mem);

    Conic::set_work(mem, arg, res, iw, w);

    m->d.prob = &p_;
    m->d.qp = &m->d_qp;

    casadi_highs_init(&m->d, &arg, &res, &iw, &w);

    for (auto&& op : opts_) {
      auto it = std::find(param_bool.begin(), param_bool.end(), op.first);
      if (it != param_bool.end() || op.second.getType() == OT_BOOL) {
        casadi_assert(kHighsStatusOk ==
          Highs_setBoolOptionValue(m->d.highs, op.first.c_str(), op.second.to_bool()),
          "Error setting option '" + op.first + "'.");
      } else if (op.second.getType() == OT_INT){
        casadi_assert(kHighsStatusOk ==
          Highs_setIntOptionValue(m->d.highs, op.first.c_str(), op.second.to_int()),
          "Error setting option '" + op.first + "'.");
      } else if (op.second.getType() == OT_DOUBLE) {
        casadi_assert(kHighsStatusOk ==
          Highs_setDoubleOptionValue(m->d.highs, op.first.c_str(), op.second.to_double()),
          "Error setting option '" + op.first + "'.");
      } else if (op.second.getType() == OT_STRING) {
        std::string v = op.second.to_string();
        casadi_assert(kHighsStatusOk ==
          Highs_setStringOptionValue(m->d.highs, op.first.c_str(), v.c_str()),
          "Error setting option '" + op.first + "'.");
      } else {
        casadi_assert(false, "Option type for '" + op.first + "'not supported!");
      }
    }

  }

  int HighsInterface::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<HighsMemory*>(mem);

    // Statistics
    m->fstats.at("solver").tic();

    casadi_highs_solve(&m->d, arg, res, iw, w);
    m->fstats.at("solver").toc();

    return 0;
  }

  HighsInterface::~HighsInterface() {
    clear_mem();
  }

  void HighsInterface::codegen_body(CodeGenerator& g) const {
    qp_codegen_body(g);
    g.add_auxiliary(CodeGenerator::AUX_PROJECT);
    g.add_auxiliary(CodeGenerator::AUX_SCAL);
    g.add_auxiliary(CodeGenerator::AUX_SPARSIFY);
    g.add_auxiliary(CodeGenerator::AUX_MAX);
    g.add_auxiliary(CodeGenerator::AUX_SPARSITY);
    g.add_auxiliary(CodeGenerator::AUX_SUM);
    g.add_auxiliary(CodeGenerator::AUX_FILL);
    g.add_auxiliary(CodeGenerator::AUX_CLIP_MIN);
    g.add_auxiliary(CodeGenerator::AUX_CLIP_MAX);
    g.add_auxiliary(CodeGenerator::AUX_DOT);
    g.add_auxiliary(CodeGenerator::AUX_BILIN);
    g.add_include("interfaces/highs_c_api.h");

    g.auxiliaries << g.sanitize_source(highs_runtime_str, {"casadi_real"});


    g.local("d", "struct casadi_highs_data*");
    g.init_local("d", "&" + codegen_mem(g));
    g.local("p", "struct casadi_highs_prob");
    set_highs_prob(g);

    // Setup data structure (corresponds to set_work)
    g << "d->prob = &p;\n";
    g << "d->qp = &d_qp;\n";
    g << "casadi_highs_init(d, &arg, &res, &iw, &w);\n";
  
    for (auto&& op : opts_) {
      auto it = std::find(param_bool.begin(), param_bool.end(), op.first);
      if (it != param_bool.end() || op.second.getType() == OT_BOOL) {
        g << "Highs_setBoolOptionValue(d->highs, " << g.constant(op.first) << ", "
          << int(op.second.to_bool()) << ");\n";
      } else if (op.second.getType() == OT_INT){
        g << "Highs_setIntOptionValue(d->highs, " << g.constant(op.first) << ", "
          << int(op.second.to_int()) << ");\n";
      } else if (op.second.getType() == OT_DOUBLE) {
        g << "Highs_setDoubleOptionValue(d->highs, " << g.constant(op.first) << ", "
          << g.constant(op.second.to_int()) << ");\n";
      } else if (op.second.getType() == OT_STRING) {
        g << "Highs_setStringOptionValue(d->highs, " << g.constant(op.first) << ", "
          << g.constant(op.second.to_string()) << ");\n";
      } else {
        casadi_assert(false, "Option type for '" + op.first + "'not supported!");
      }
    }

    g << "casadi_highs_solve(d, arg, res, iw, w);\n";

    g << "if (!d_qp.success) {\n";
    if (error_on_fail_) {
      g << "return -1000;\n";
    } else {
      g << "return -1;\n";
    }
    g << "}\n";
    g << "return 0;\n";
  }

  Dict HighsInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    Highs highs;
    auto m = static_cast<HighsMemory*>(mem);
    stats["return_status"] = highs.modelStatusToString(static_cast<HighsModelStatus>(m->d.return_status));
    stats["simplex_iteration_count"] = m->d.simplex_iteration_count;
    stats["simplex_iteration_count"] = m->d.simplex_iteration_count;
    stats["ipm_iteration_count"] = m->d.ipm_iteration_count;
    stats["qp_iteration_count"] = m->d.qp_iteration_count;
    stats["crossover_iteration_count"] = m->d.crossover_iteration_count;
    stats["primal_solution_status"]  = highs.solutionStatusToString(m->d.primal_solution_status);
    stats["dual_solution_status"] = highs.solutionStatusToString(m->d.dual_solution_status);
    stats["basis_validity"] = highs.basisValidityToString(m->d.basis_validity);
    stats["mip_dual_bound"] = m->d.mip_dual_bound;
    stats["mip_gap"] = m->d.mip_gap;
    stats["num_primal_infeasibilities"] = m->d.num_primal_infeasibilities;
    stats["max_primal_infeasibility"] = m->d.max_primal_infeasibility;
    stats["sum_primal_infeasibilities"] = m->d.sum_primal_infeasibilities;
    stats["num_dual_infeasibilities"] = m->d.num_dual_infeasibilities;
    stats["max_dual_infeasibility"] = m->d.max_dual_infeasibility;
    stats["sum_dual_infeasibilities"] = m->d.sum_dual_infeasibilities;
    return stats;
  }

  HighsInterface::HighsInterface(DeserializingStream& s) : Conic(s) {
    s.version("HighsInterface", 1);
    s.unpack("HighsInterface::opts", opts_);
    init_dependent();
    set_highs_prob();
  }

  void HighsInterface::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);

    s.version("HighsInterface", 1);
    s.pack("HighsInterface::opts", opts_);
  }

} // end namespace casadi
