/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2026 Joel Andersson, Joris Gillis, Moritz Diehl,
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

#include "mosek_interface.hpp"
#include "casadi/core/nlp_tools.hpp"
#include "casadi/core/casadi_misc.hpp"

#include <mosek_runtime_str.h>

#include <cstdio>
#include <cstring>

namespace casadi {

  extern "C"
  int CASADI_CONIC_MOSEK_EXPORT
  casadi_register_conic_mosek(Conic::Plugin* plugin) {
    plugin->creator = MosekInterface::creator;
    plugin->name = "mosek";
    plugin->doc = MosekInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &MosekInterface::options_;
    plugin->deserialize = &MosekInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_MOSEK_EXPORT casadi_load_conic_mosek() {
    Conic::registerPlugin(casadi_register_conic_mosek);
  }


  MosekInterface::MosekInterface(const std::string& name,
                                 const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {
  }

  const Options MosekInterface::options_
  = {{&Conic::options_},
     {{"mosek",
       {OT_DICT,
        "Options to be passed to MOSEK.  Each entry's key is the MOSEK "
        "parameter name (e.g. \"MSK_IPAR_LOG\", \"MSK_DPAR_INTPNT_QO_TOL_REL_GAP\"); "
        "the value's type must match the parameter's type "
        "(int parameters take int, dpar take double, spar take string)."
        }},
     }
   };

  void MosekInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Conic::init(opts);

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="mosek") {
        opts_ = op.second;
      }
    }

    // Pre-compute SOCP mapping (must run before set_mosek_prob so the
    // socp pointer is resolved in p_)
    has_socp_ = !Q_.is_null() && Q_.nnz() > 0;
    if (has_socp_) build_socp_config();

    init_dependent();
    set_mosek_prob();

    // Allocate memory
    casadi_int sz_arg, sz_res, sz_w, sz_iw;
    casadi_mosek_work(&p_, &sz_arg, &sz_res, &sz_iw, &sz_w);

    alloc_arg(sz_arg, true);
    alloc_res(sz_res, true);
    alloc_iw(sz_iw, true);
    alloc_w(sz_w, true);
  }

  void MosekInterface::init_dependent() {
    colinda_.resize(A_.size2() + 1);
    rowa_.resize(A_.nnz());
    copy_vector(A_.colind(), colinda_);
    copy_vector(A_.row(), rowa_);

    // Build lower-triangular triplets from the symmetric H_.
    // Mosek's MSK_putqobj wants entries with row >= col (lower triangle).
    qobj_row_.clear();
    qobj_col_.clear();
    qobj_nz_idx_.clear();
    const casadi_int* Hc = H_.colind();
    const casadi_int* Hr = H_.row();
    for (casadi_int j = 0; j < H_.size2(); ++j) {
      for (casadi_int k = Hc[j]; k < Hc[j + 1]; ++k) {
        casadi_int i = Hr[k];
        if (i >= j) {
          qobj_row_.push_back(static_cast<int>(i));
          qobj_col_.push_back(static_cast<int>(j));
          qobj_nz_idx_.push_back(static_cast<int>(k));
        }
      }
    }

    // Discrete variable flags: 'I' for integer, 'C' for continuous.
    // Only mark as MIP (non-empty coltype_) if at least one variable is
    // actually discrete -- callers (like Opti) pass an all-false discrete
    // vector for pure LPs/QPs, which we should NOT classify as MIP.
    coltype_.clear();
    bool any_discrete = false;
    for (bool d : discrete_) if (d) { any_discrete = true; break; }
    if (any_discrete) {
      coltype_.assign(nx_, 'C');
      for (casadi_int i = 0; i < nx_; ++i) {
        if (discrete_[i]) coltype_[i] = 'I';
      }
    }
  }

  void MosekInterface::set_mosek_prob() {
    p_.qp = &p_qp_;
    p_.colinda = get_ptr(colinda_);
    p_.rowa = get_ptr(rowa_);
    p_.qobj_row = qobj_row_.empty() ? nullptr : get_ptr(qobj_row_);
    p_.qobj_col = qobj_col_.empty() ? nullptr : get_ptr(qobj_col_);
    p_.qobj_nz_idx = qobj_nz_idx_.empty() ? nullptr : get_ptr(qobj_nz_idx_);
    p_.nquad = static_cast<int>(qobj_row_.size());
    p_.coltype = coltype_.empty() ? nullptr : get_ptr(coltype_);

    p_.socp = has_socp_ ? &socp_ : nullptr;

    casadi_mosek_setup(&p_);
  }

  // Emit "p.<field> = <local int[]>;" for an int*-typed prob field.
  // Empty vectors emit a NULL assignment directly (avoids the unused-constant
  // warning from constant_copy on empty arrays).
  static void codegen_int_field(CodeGenerator& g, const std::string& prob_field,
      const std::string& local_name, const std::vector<int>& v) {
    if (v.empty()) {
      g << "p." << prob_field << " = 0;\n";
    } else {
      g.constant_copy(local_name, vector_static_cast<casadi_int>(v), "int");
      g << "p." << prob_field << " = " << local_name << ";\n";
    }
  }

  void MosekInterface::set_mosek_prob(CodeGenerator& g) const {
    g << "p.qp = &p_qp;\n";

    codegen_int_field(g, "colinda",     "colinda",     colinda_);
    codegen_int_field(g, "rowa",        "rowa",        rowa_);
    codegen_int_field(g, "qobj_row",    "qobj_row",    qobj_row_);
    codegen_int_field(g, "qobj_col",    "qobj_col",    qobj_col_);
    codegen_int_field(g, "qobj_nz_idx", "qobj_nz_idx", qobj_nz_idx_);
    g << "p.nquad = " << qobj_row_.size() << ";\n";

    if (coltype_.empty()) {
      g << "p.coltype = 0;\n";
    } else {
      g << "p.coltype = " << g.constant(coltype_) << ";\n";
    }

    if (has_socp_) {
      g.local("socp_p", "struct casadi_socp_prob");
      g << "socp_p.n_blocks = " << socp_.n_blocks << ";\n";
      g << "socp_p.n_lifted = " << socp_.n_lifted << ";\n";
      g << "socp_p.nx = " << socp_.nx << ";\n";
      g << "socp_p.r = " << g.constant(socp_r_) << ";\n";
      g << "socp_p.mq_colind = " << g.constant(socp_mq_colind_) << ";\n";
      g << "socp_p.mq_row = " << g.constant(socp_mq_row_) << ";\n";
      g << "socp_p.mq_data = " << g.constant(socp_mq_data_) << ";\n";
      g << "socp_p.map_P = " << g.constant(socp_map_P_) << ";\n";
      g << "casadi_socp_setup(&socp_p);\n";
      g << "p.socp = &socp_p;\n";
    } else {
      g << "p.socp = 0;\n";
    }

    g << "casadi_mosek_setup(&p);\n";
  }

  void MosekInterface::codegen_init_mem(CodeGenerator& g) const {
    g << "casadi_mosek_init_mem(&" + codegen_mem(g) + ");\n";
    g << "return 0;\n";
  }

  void MosekInterface::codegen_free_mem(CodeGenerator& g) const {
    g << "casadi_mosek_free_mem(&" + codegen_mem(g) + ");\n";
  }

  void MosekInterface::codegen_body(CodeGenerator& g) const {
    qp_codegen_body(g);
    g.add_auxiliary(CodeGenerator::AUX_INF);
    // Always emit casadi_socp definitions: the mosek data struct
    // references casadi_socp_data even when no Q/P is present.
    g.add_auxiliary(CodeGenerator::AUX_SOCP);
    g.add_include("mosek.h");
    g.auxiliaries << g.sanitize_source(mosek_runtime_str, {"casadi_real"});

    g.local("d", "struct casadi_mosek_data*");
    g.init_local("d", "&" + codegen_mem(g));
    g.local("p", "struct casadi_mosek_prob");
    set_mosek_prob(g);

    g << "d->prob = &p;\n";
    g << "d->qp = &d_qp;\n";
    g << "casadi_mosek_init(d, &arg, &res, &iw, &w);\n";

    if (has_socp_) {
      g << "d->socp.q = arg[" << CONIC_Q << "];\n";
      g << "d->socp.p = arg[" << CONIC_P << "];\n";
    }

    // Apply user options.  Resolve the parameter's type at runtime via
    // MSK_getparamname/MSK_getparamtype is awkward; the simplest robust
    // path is to dispatch on the CasADi-side value type and let Mosek's
    // typed setter accept-or-reject by name.
    for (auto&& op : opts_) {
      if (op.second.is_double()) {
        g << "MSK_putnadouparam(d->task, " << g.constant(op.first) << ", "
          << op.second.to_double() << ");\n";
      } else if (op.second.is_int() || op.second.is_bool()) {
        g << "MSK_putnaintparam(d->task, " << g.constant(op.first) << ", "
          << static_cast<int>(op.second.to_int()) << ");\n";
      } else if (op.second.is_string()) {
        g << "MSK_putnastrparam(d->task, " << g.constant(op.first) << ", "
          << g.constant(op.second.to_string()) << ");\n";
      } else {
        casadi_error("Unsupported option type for '" + op.first + "'.");
      }
    }

    g << "casadi_mosek_solve(d, arg, res, iw, w);\n";

    g << "if (!d_qp.success) {\n";
    if (error_on_fail_) {
      g << "  return -1000;\n";
    } else {
      g << "  return -1;\n";
    }
    g << "}\n";
    g << "return 0;\n";
  }

  void MosekInterface::build_socp_config() {
    SDPToSOCPMem sm;
    sdp_to_socp_init(sm);

    socp_r_ = sm.r;
    casadi_int n_eq = sm.map_Q.size2();
    socp_mq_colind_.assign(sm.map_Q.colind(), sm.map_Q.colind() + n_eq + 1);
    socp_mq_row_.assign(sm.map_Q.row(), sm.map_Q.row() + sm.map_Q.nnz());
    socp_mq_data_.assign(sm.map_Q.ptr(), sm.map_Q.ptr() + sm.map_Q.nnz());
    socp_map_P_ = sm.map_P;

    socp_.n_blocks = static_cast<casadi_int>(sm.r.size() - 1);
    socp_.r = get_ptr(socp_r_);
    socp_.n_lifted = sm.r.back();
    socp_.nx = nx_;
    socp_.mq_colind = get_ptr(socp_mq_colind_);
    socp_.mq_row = get_ptr(socp_mq_row_);
    socp_.mq_data = get_ptr(socp_mq_data_);
    socp_.map_P = get_ptr(socp_map_P_);
    casadi_socp_setup(&socp_);
  }

  static void mosek_stream_cb(void* /*h*/, const char* s) { std::fputs(s, stderr); }

  int MosekInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    if (!mem) return 1;

    auto m = static_cast<MosekMemory*>(mem);
    if (casadi_mosek_init_mem(&m->d)) {
      casadi_error("MOSEK: MSK_makeenv()/MSK_maketask() failed. "
        "Check that the MOSEK license is reachable "
        "(MOSEKLM_LICENSE_FILE or $HOME/mosek/mosek.lic).");
    }
    if (verbose_) {
      MSK_linkfunctotaskstream(m->d.task, MSK_STREAM_LOG, nullptr, mosek_stream_cb);
    }

    m->add_stat("preprocessing");
    m->add_stat("solver");
    m->add_stat("postprocessing");

    return 0;
  }

  void MosekInterface::free_mem(void* mem) const {
    auto m = static_cast<MosekMemory*>(mem);
    casadi_mosek_free_mem(&m->d);
    delete static_cast<MosekMemory*>(mem);
  }

  /** \brief Set the (persistent) work vectors */
  void MosekInterface::set_work(void* mem, const double**& arg, double**& res,
                                casadi_int*& iw, double*& w) const {

    auto m = static_cast<MosekMemory*>(mem);

    Conic::set_work(mem, arg, res, iw, w);

    m->d.prob = &p_;
    m->d.qp = &m->d_qp;

    casadi_mosek_init(&m->d, &arg, &res, &iw, &w);

    // SOCP inputs
    if (has_socp_) {
      m->d.socp.q = arg[CONIC_Q];
      m->d.socp.p = arg[CONIC_P];
    }

    // Push user parameters to MOSEK using typed setters dispatched on
    // CasADi-side value type.  The "na" (name-addressed) variants accept
    // a name string and silently fail with a return code if the name does
    // not match the type -- assert success here so the user gets a clear
    // error.
    for (auto&& op : opts_) {
      MSKrescodee rc = MSK_RES_OK;
      if (op.second.is_double()) {
        rc = MSK_putnadouparam(m->d.task, op.first.c_str(), op.second.to_double());
      } else if (op.second.is_int() || op.second.is_bool()) {
        rc = MSK_putnaintparam(m->d.task, op.first.c_str(),
                               static_cast<int>(op.second.to_int()));
      } else if (op.second.is_string()) {
        rc = MSK_putnastrparam(m->d.task, op.first.c_str(),
                               op.second.to_string().c_str());
      } else {
        casadi_error("MOSEK: unsupported value type for parameter '"
                     + op.first + "'.");
      }
      casadi_assert(rc == MSK_RES_OK,
        "MOSEK: failed to set parameter '" + op.first + "' (rc=" + str(rc) + ").");
    }
  }

  int MosekInterface::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<MosekMemory*>(mem);

    m->fstats.at("solver").tic();
    int rc = casadi_mosek_solve(&m->d, arg, res, iw, w);
    m->fstats.at("solver").toc();

    return rc;
  }

  MosekInterface::~MosekInterface() {
    clear_mem();
  }

  Dict MosekInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<MosekMemory*>(mem);
    stats["return_status"] = m->d.return_status;
    stats["prob_status"] = m->d.prob_status;
    stats["sol_status"] = m->d.sol_status;
    stats["objective"] = m->d.obj_val;
    return stats;
  }

  MosekInterface::MosekInterface(DeserializingStream& s) : Conic(s) {
    s.version("MosekInterface", 1);
    s.unpack("MosekInterface::opts", opts_);
    has_socp_ = !Q_.is_null() && Q_.nnz() > 0;
    if (has_socp_) build_socp_config();
    init_dependent();
    set_mosek_prob();
  }

  void MosekInterface::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);

    s.version("MosekInterface", 1);
    s.pack("MosekInterface::opts", opts_);
  }

} // end namespace casadi
