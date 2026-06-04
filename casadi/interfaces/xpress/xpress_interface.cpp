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

#include "xpress_interface.hpp"
#include "casadi/core/nlp_tools.hpp"
#include "casadi/core/casadi_misc.hpp"

#include <xpress_runtime_str.h>

#include <mutex>

namespace casadi {

  // Xpress requires a single global library init via XPRSinit().
  // We do it lazily on first plugin use; XPRSfree() is registered atexit
  // by the OS-level library teardown when the process ends.
  static std::once_flag s_xpress_init_flag;
  static int s_xpress_init_rc = 0;
  static void xpress_global_init() {
    s_xpress_init_rc = XPRSinit(nullptr);
  }

  extern "C"
  int CASADI_CONIC_XPRESS_EXPORT
  casadi_register_conic_xpress(Conic::Plugin* plugin) {
    plugin->creator = XpressInterface::creator;
    plugin->name = "xpress";
    plugin->doc = XpressInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &XpressInterface::options_;
    plugin->deserialize = &XpressInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_XPRESS_EXPORT casadi_load_conic_xpress() {
    Conic::registerPlugin(casadi_register_conic_xpress);
  }


  XpressInterface::XpressInterface(const std::string& name,
                                   const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {
  }

  const Options XpressInterface::options_
  = {{&Conic::options_},
     {{"xpress",
       {OT_DICT,
        "Options to be passed to FICO Xpress.  Each entry's key is the "
        "Xpress control name (e.g. \"OUTPUTLOG\", \"MAXTIME\", \"MIPRELSTOP\"); "
        "the value's type must match the control's type (int / double / string)."
        }},
      {"sos_groups",
       {OT_INTVECTORVECTOR,
        "Definition of SOS groups by indices."}},
      {"sos_weights",
       {OT_DOUBLEVECTORVECTOR,
        "Weights corresponding to SOS entries."}},
      {"sos_types",
       {OT_INTVECTOR,
        "Specify 1 or 2 for each SOS group."}},
     }
   };

  void XpressInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Conic::init(opts);

    std::vector< std::vector<casadi_int> > sos_groups;
    std::vector< std::vector<double> > sos_weights;
    std::vector<casadi_int> sos_types;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="xpress") {
        opts_ = op.second;
      } else if (op.first=="sos_groups") {
        sos_groups = op.second.to_int_vector_vector();
      } else if (op.first=="sos_weights") {
        sos_weights = op.second.to_double_vector_vector();
      } else if (op.first=="sos_types") {
        sos_types = op.second.to_int_vector();
      }
    }

    // Validate + populate defaults (weights default to 1..n, types default to 1)
    check_sos(nx_, sos_groups, sos_weights, sos_types);

    // Flatten into Xpress's CSR-like storage
    if (!sos_groups.empty()) {
      std::vector<casadi_int> beg, ind;
      flatten_nested_vector(sos_groups, ind, beg);
      flatten_nested_vector(sos_weights, sos_refval_);
      sos_setstart_.assign(beg.begin(), beg.end());
      sos_setind_.assign(ind.begin(), ind.end());
      sos_settype_.resize(sos_types.size());
      for (size_t k = 0; k < sos_types.size(); ++k) {
        sos_settype_[k] = static_cast<char>('0' + sos_types[k]);  // '1' or '2'
      }
    }

    // Pre-compute SOCP mapping (must run before set_xpress_prob so the
    // socp pointer is resolved in p_)
    has_socp_ = !Q_.is_null() && Q_.nnz() > 0;
    if (has_socp_) build_socp_config();

    init_dependent();
    set_xpress_prob();

    // Allocate memory
    casadi_int sz_arg, sz_res, sz_w, sz_iw;
    casadi_xpress_work(&p_, &sz_arg, &sz_res, &sz_iw, &sz_w);

    alloc_arg(sz_arg, true);
    alloc_res(sz_res, true);
    alloc_iw(sz_iw, true);
    alloc_w(sz_w, true);
  }

  void XpressInterface::init_dependent() {
    colinda_.resize(A_.size2() + 1);
    rowa_.resize(A_.nnz());
    copy_vector(A_.colind(), colinda_);
    copy_vector(A_.row(), rowa_);

    // Build upper-triangular triplets from the symmetric H_.
    // H_ is stored CSC; iterate columns and pick (row, col) with row <= col.
    qobj_col1_.clear();
    qobj_col2_.clear();
    qobj_nz_idx_.clear();
    const casadi_int* Hc = H_.colind();
    const casadi_int* Hr = H_.row();
    for (casadi_int j = 0; j < H_.size2(); ++j) {
      for (casadi_int k = Hc[j]; k < Hc[j + 1]; ++k) {
        casadi_int i = Hr[k];
        if (i <= j) {
          qobj_col1_.push_back(static_cast<int>(i));
          qobj_col2_.push_back(static_cast<int>(j));
          qobj_nz_idx_.push_back(static_cast<int>(k));
        }
      }
    }

    // Discrete variable flags: 'I' for integer, 'C' for continuous
    if (!discrete_.empty()) {
      coltype_.assign(nx_, 'C');
      for (casadi_int i = 0; i < nx_; ++i) {
        if (discrete_[i]) coltype_[i] = 'I';
      }
    }
  }

  void XpressInterface::set_xpress_prob() {
    p_.qp = &p_qp_;
    p_.colinda = get_ptr(colinda_);
    p_.rowa = get_ptr(rowa_);
    p_.qobj_col1 = get_ptr(qobj_col1_);
    p_.qobj_col2 = get_ptr(qobj_col2_);
    p_.qobj_nz_idx = get_ptr(qobj_nz_idx_);
    p_.nquad = static_cast<int>(qobj_col1_.size());
    p_.coltype = coltype_.empty() ? nullptr : get_ptr(coltype_);

    p_.n_sos_sets = static_cast<int>(sos_settype_.size());
    p_.n_sos_elems = static_cast<int>(sos_setind_.size());
    p_.sos_settype = sos_settype_.empty() ? nullptr : get_ptr(sos_settype_);
    p_.sos_setstart = sos_setstart_.empty() ? nullptr : get_ptr(sos_setstart_);
    p_.sos_setind = sos_setind_.empty() ? nullptr : get_ptr(sos_setind_);
    p_.sos_refval = sos_refval_.empty() ? nullptr : get_ptr(sos_refval_);

    p_.socp = has_socp_ ? &socp_ : nullptr;

    casadi_xpress_setup(&p_);
  }

  // Emit "p.<field> = <local int[]>;" for an int*-typed prob field.
  // For empty vectors, emit a NULL assignment directly instead of going
  // through constant_copy (which emits an unused constant ref).
  static void codegen_int_field(CodeGenerator& g, const std::string& prob_field,
      const std::string& local_name, const std::vector<int>& v) {
    if (v.empty()) {
      g << "p." << prob_field << " = 0;\n";
    } else {
      g.constant_copy(local_name, vector_static_cast<casadi_int>(v), "int");
      g << "p." << prob_field << " = " << local_name << ";\n";
    }
  }

  void XpressInterface::set_xpress_prob(CodeGenerator& g) const {
    g << "p.qp = &p_qp;\n";

    codegen_int_field(g, "colinda",     "colinda",     colinda_);
    codegen_int_field(g, "rowa",        "rowa",        rowa_);
    codegen_int_field(g, "qobj_col1",   "qobj_col1",   qobj_col1_);
    codegen_int_field(g, "qobj_col2",   "qobj_col2",   qobj_col2_);
    codegen_int_field(g, "qobj_nz_idx", "qobj_nz_idx", qobj_nz_idx_);
    g << "p.nquad = " << qobj_col1_.size() << ";\n";

    if (coltype_.empty()) {
      g << "p.coltype = 0;\n";
    } else {
      g << "p.coltype = " << g.constant(coltype_) << ";\n";
    }

    if (sos_settype_.empty()) {
      g << "p.n_sos_sets = 0;\n";
      g << "p.n_sos_elems = 0;\n";
      g << "p.sos_settype = 0;\n";
      g << "p.sos_setstart = 0;\n";
      g << "p.sos_setind = 0;\n";
      g << "p.sos_refval = 0;\n";
    } else {
      codegen_int_field(g, "sos_setstart", "sos_setstart", sos_setstart_);
      codegen_int_field(g, "sos_setind",   "sos_setind",   sos_setind_);
      g << "p.sos_settype = " << g.constant(sos_settype_) << ";\n";
      g << "p.sos_refval = " << g.constant(sos_refval_) << ";\n";
      g << "p.n_sos_sets = " << sos_settype_.size() << ";\n";
      g << "p.n_sos_elems = " << sos_setind_.size() << ";\n";
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

    g << "casadi_xpress_setup(&p);\n";
  }

  void XpressInterface::codegen_init_mem(CodeGenerator& g) const {
    g << "casadi_xpress_init_mem(&" + codegen_mem(g) + ");\n";
    g << "return 0;\n";
  }

  void XpressInterface::codegen_free_mem(CodeGenerator& g) const {
    g << "casadi_xpress_free_mem(&" + codegen_mem(g) + ");\n";
  }

  void XpressInterface::codegen_body(CodeGenerator& g) const {
    qp_codegen_body(g);
    g.add_auxiliary(CodeGenerator::AUX_INF);
    // Always emit casadi_socp definitions: the xpress data struct
    // references casadi_socp_data even when no Q/P is present.
    g.add_auxiliary(CodeGenerator::AUX_SOCP);
    g.add_include("xprs.h");
    g.auxiliaries << g.sanitize_source(xpress_runtime_str, {"casadi_real"});

    g.local("d", "struct casadi_xpress_data*");
    g.init_local("d", "&" + codegen_mem(g));
    g.local("p", "struct casadi_xpress_prob");
    set_xpress_prob(g);

    g << "d->prob = &p;\n";
    g << "d->qp = &d_qp;\n";
    g << "casadi_xpress_init(d, &arg, &res, &iw, &w);\n";

    if (has_socp_) {
      g << "d->socp.q = arg[" << CONIC_Q << "];\n";
      g << "d->socp.p = arg[" << CONIC_P << "];\n";
    }

    // Set crossover default (matches C++ set_work)
    g << "XPRSsetintcontrol(d->xprob, " << XPRS_CROSSOVER << ", 1);\n";

    // Apply user options.  Dispatch by the CasADi-side type tag of each
    // value; the control id is resolved at runtime via XPRSgetcontrolinfo
    // (which lets the generated code work without an Xpress license at
    // codegen time).
    for (auto&& op : opts_) {
      g << "{\n";
      g << "  int xprs_id, xprs_type;\n";
      g << "  XPRSgetcontrolinfo(d->xprob, " << g.constant(op.first)
        << ", &xprs_id, &xprs_type);\n";
      if (op.second.is_double()) {
        g << "  XPRSsetdblcontrol(d->xprob, xprs_id, "
          << op.second.to_double() << ");\n";
      } else if (op.second.is_int() || op.second.is_bool()) {
        g << "  XPRSsetintcontrol(d->xprob, xprs_id, "
          << static_cast<int>(op.second.to_int()) << ");\n";
      } else if (op.second.is_string()) {
        g << "  XPRSsetstrcontrol(d->xprob, xprs_id, "
          << g.constant(op.second.to_string()) << ");\n";
      } else {
        casadi_error("Unsupported option type for '" + op.first + "'.");
      }
      g << "}\n";
    }

    g << "casadi_xpress_solve(d, arg, res, iw, w);\n";

    g << "if (!d_qp.success) {\n";
    if (error_on_fail_) {
      g << "  return -1000;\n";
    } else {
      g << "  return -1;\n";
    }
    g << "}\n";
    g << "return 0;\n";
  }

  void XpressInterface::build_socp_config() {
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

  int XpressInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    if (!mem) return 1;

    // Lazy global library init
    std::call_once(s_xpress_init_flag, xpress_global_init);
    casadi_assert(s_xpress_init_rc == 0,
      "XPRSinit() failed (return code " + str(s_xpress_init_rc) + "). "
      "Check that the XPRESS environment variable points to a valid license "
      "or that a license server is reachable.");

    auto m = static_cast<XpressMemory*>(mem);
    if (casadi_xpress_init_mem(&m->d)) return 1;

    m->add_stat("preprocessing");
    m->add_stat("solver");
    m->add_stat("postprocessing");

    return 0;
  }

  void XpressInterface::free_mem(void* mem) const {
    auto m = static_cast<XpressMemory*>(mem);
    casadi_xpress_free_mem(&m->d);
    delete static_cast<XpressMemory*>(mem);
  }

  /** \brief Set the (persistent) work vectors */
  void XpressInterface::set_work(void* mem, const double**& arg, double**& res,
                                 casadi_int*& iw, double*& w) const {

    auto m = static_cast<XpressMemory*>(mem);

    Conic::set_work(mem, arg, res, iw, w);

    m->d.prob = &p_;
    m->d.qp = &m->d_qp;

    casadi_xpress_init(&m->d, &arg, &res, &iw, &w);

    // SOCP inputs
    if (has_socp_) {
      m->d.socp.q = arg[CONIC_Q];
      m->d.socp.p = arg[CONIC_P];
    }

    // Sensible default: enable crossover after the QP barrier so the
    // returned solution sits exactly at the optimal vertex.  Without this
    // Xpress returns an interior-point near-vertex solution that can be
    // off by ~FEASTOL on each active constraint.  User options below
    // can still override this.
    XPRSsetintcontrol(m->d.xprob, XPRS_CROSSOVER, 1);

    // Push user options to Xpress.  Resolve name -> (id, type) via
    // XPRSgetcontrolinfo, then dispatch to the typed setter.
    for (auto&& op : opts_) {
      int id = 0;
      int type = 0;
      int rc = XPRSgetcontrolinfo(m->d.xprob, op.first.c_str(), &id, &type);
      casadi_assert(rc == 0,
        "Xpress: unknown control '" + op.first + "' (XPRSgetcontrolinfo rc=" +
        str(rc) + ").");
      switch (type) {
        case XPRS_TYPE_INT:
          rc = XPRSsetintcontrol(m->d.xprob, id, op.second.to_int());
          break;
        case XPRS_TYPE_INT64:
          rc = XPRSsetintcontrol64(m->d.xprob, id,
              static_cast<XPRSint64>(op.second.to_int()));
          break;
        case XPRS_TYPE_DOUBLE:
          rc = XPRSsetdblcontrol(m->d.xprob, id, op.second.to_double());
          break;
        case XPRS_TYPE_STRING: {
          std::string v = op.second.to_string();
          rc = XPRSsetstrcontrol(m->d.xprob, id, v.c_str());
          break;
        }
        default:
          casadi_error("Xpress: unsupported control type for '" + op.first + "'.");
      }
      casadi_assert(rc == 0,
        "Xpress: failed to set control '" + op.first + "' (rc=" + str(rc) + ").");
    }
  }

  int XpressInterface::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<XpressMemory*>(mem);

    m->fstats.at("solver").tic();
    int rc = casadi_xpress_solve(&m->d, arg, res, iw, w);
    m->fstats.at("solver").toc();

    return rc;
  }

  XpressInterface::~XpressInterface() {
    clear_mem();
  }

  Dict XpressInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<XpressMemory*>(mem);
    stats["return_status"] = m->d.return_status;
    stats["lp_status"] = m->d.lp_status;
    stats["mip_status"] = m->d.mip_status;
    stats["simplex_iter"] = m->d.simplex_iter;
    stats["barrier_iter"] = m->d.barrier_iter;
    stats["mip_nodes"] = m->d.mip_nodes;
    stats["objective"] = m->d.obj_val;
    return stats;
  }

  XpressInterface::XpressInterface(DeserializingStream& s) : Conic(s) {
    s.version("XpressInterface", 1);
    s.unpack("XpressInterface::opts", opts_);
    s.unpack("XpressInterface::sos_settype", sos_settype_);
    s.unpack("XpressInterface::sos_setstart", sos_setstart_);
    s.unpack("XpressInterface::sos_setind", sos_setind_);
    s.unpack("XpressInterface::sos_refval", sos_refval_);
    has_socp_ = !Q_.is_null() && Q_.nnz() > 0;
    if (has_socp_) build_socp_config();
    init_dependent();
    set_xpress_prob();
  }

  void XpressInterface::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);

    s.version("XpressInterface", 1);
    s.pack("XpressInterface::opts", opts_);
    s.pack("XpressInterface::sos_settype", sos_settype_);
    s.pack("XpressInterface::sos_setstart", sos_setstart_);
    s.pack("XpressInterface::sos_setind", sos_setind_);
    s.pack("XpressInterface::sos_refval", sos_refval_);
  }

} // end namespace casadi
