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

#include "conopt_interface.hpp"
#include "casadi/core/casadi_misc.hpp"
#include "casadi/core/casadi_interrupt.hpp"
#include "casadi/core/sx_function.hpp"
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <limits>

namespace casadi {

  // CONOPT 4.39.2 accepts option names up to 20 characters.
  static const size_t conopt_max_option_name = 20;

  extern "C" int CASADI_NLPSOL_CONOPT_EXPORT casadi_register_nlpsol_conopt(Nlpsol::Plugin* plugin) {
    plugin->creator = ConoptInterface::creator;
    plugin->name = "conopt";
    plugin->doc = ConoptInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &ConoptInterface::options_;
    plugin->deserialize = &ConoptInterface::deserialize;
    return 0;
  }

  extern "C" void CASADI_NLPSOL_CONOPT_EXPORT casadi_load_nlpsol_conopt() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_conopt);
  }

  const Options ConoptInterface::options_ =
      {{&Nlpsol::options_}, {
      {"exact_hessian", {OT_BOOL, "Provide exact Hessian to CONOPT"}},
      {"warm_start", {OT_BOOL,
        "Warm-start CONOPT using multipliers from a prior solve to infer "
        "basis status (IniStat=2)"}},
      {"conopt", {OT_DICT, "Options to be passed to CONOPT"}},
      {"optfile", {OT_STRING,
        "Path to a CONOPT option file (for string-valued CR-cells such as Algorithm)"}},
      {"debug", {OT_BOOL,
        "Print debug output: constraint values at each FDEval, solution vector, "
        "and option echo"}},
      {"license", {OT_DICT,
        "CONOPT license as a dict with keys 'int_1', 'int_2', 'int_3' (int) and "
        "'text' (string). Overrides the CONOPT_LICENSE_* environment variables. "
        "Not serialized."}},
      {"subset_eval", {OT_BOOL,
        "Evaluate only the rows CONOPT lists in FDEvalIni, by executing just the "
        "instructions those rows depend on. Requires the evaluation functions to be "
        "SX without function calls (e.g. set 'expand' to true for MX problems); "
        "otherwise all rows are evaluated. Default: true"}}
  }};

  ConoptInterface::ConoptInterface(const std::string& name, const Function& nlp)
      : Nlpsol(name, nlp) {}

  ConoptInterface::~ConoptInterface() { clear_mem(); }

  void ConoptInterface::init(const Dict& opts) {
    Nlpsol::init(opts);

    // Extract native options
    warm_start_ = false;
    debug_ = false;
    subset_eval_ = true;
    for (auto&& op : opts) {
      if (op.first == "conopt") opts_ = op.second;
      else if (op.first == "optfile") optfile_ = op.second.to_string();
      else if (op.first == "warm_start") warm_start_ = op.second.to_bool();
      else if (op.first == "debug") debug_ = op.second.to_bool();
      else if (op.first == "subset_eval") subset_eval_ = op.second.to_bool();
      else if (op.first == "license") {
        Dict lic = op.second;
        for (const char* key : {"int_1", "int_2", "int_3", "text"}) {
          casadi_assert(lic.find(key) != lic.end(),
            "CONOPT 'license' option is missing key '" + std::string(key) + "'. "
            "Expected keys: int_1, int_2, int_3, text.");
        }
        for (auto&& e : lic) {
          if (e.first == "int_1") license_int_[0] = e.second.to_int();
          else if (e.first == "int_2") license_int_[1] = e.second.to_int();
          else if (e.first == "int_3") license_int_[2] = e.second.to_int();
          else if (e.first == "text") license_text_ = e.second.to_string();
          else casadi_error("CONOPT 'license' option: unknown key '" + e.first + "'.");
        }
        has_license_ = true;
      }
    }

    for (auto&& op : opts_) {
      casadi_assert(op.first.size() <= conopt_max_option_name,
        "CONOPT option name '" + op.first + "' is " + str(op.first.size()) +
        " characters; CONOPT accepts at most " + str(conopt_max_option_name) + ".");
    }

    Function gradf_fcn = create_function("nlp_grad_f", {"x", "p"}, {"f", "grad:f:x"});
    gradf_sp_ = gradf_fcn.sparsity_out(1);

    Function jacg_fcn = create_function("nlp_jac_g", {"x", "p"}, {"g", "jac:g:x"});
    jacg_sp_ = jacg_fcn.sparsity_out(1);

    // Value-only function for MODE=1 FDEvalIni calls, which CONOPT issues often
    // (line searches, preprocessing); cheaper than the derivative functions above.
    // f and g together, so shared subexpressions (e.g. integrator calls in
    // shooting formulations) are evaluated once.
    create_function("nlp_fg", {"x", "p"}, {"f", "g"});

    // Same for derivative FDEvalIni calls: values and first derivatives in one
    // function. nlp_grad_f/nlp_jac_g are kept for the linearity detection and the
    // one-off evaluations at x0 in solve().
    Function fg_jac_fcn = create_function("nlp_fg_jac", {"x", "p"},
                                          {"f", "grad:f:x", "g", "jac:g:x"});
    casadi_assert(fg_jac_fcn.sparsity_out(1) == gradf_sp_ &&
                  fg_jac_fcn.sparsity_out(3) == jacg_sp_,
                  "nlp_fg_jac sparsity differs from nlp_grad_f/nlp_jac_g");

    // Detect linear (constant) Jacobian entries using second-order sparsity:
    // d(jac:g:x_compact)/dx has shape (nnz_g, nx_); if row k is empty,
    // the k-th Jacobian nonzero is constant in x (linear entry).
    {
      Sparsity djac_dx = jacg_fcn.sparsity_jac(0, 1, true);
      jacg_nlflag_.assign(jacg_sp_.nnz(), 0);
      const casadi_int* dj_row = djac_dx.row();
      for (casadi_int el = 0; el < djac_dx.nnz(); ++el)
        jacg_nlflag_[dj_row[el]] = 1;
      has_linear_jac_ = std::any_of(jacg_nlflag_.begin(), jacg_nlflag_.end(),
                                     [](int f) { return f == 0; });
    }

    // Setup 2nd Order Info
    exact_hessian_ = true;
    if (opts.find("exact_hessian") != opts.end()) exact_hessian_ = opts.at("exact_hessian");

    if (exact_hessian_) {
      Function hl_fcn = create_function("nlp_hess_l", {"x", "p", "lam:f", "lam:g"},
                                        {"tril:hess:gamma:x:x"}, {{"gamma", {"f", "g"}}});
      hesslag_sp_ = hl_fcn.sparsity_out(0);
    }

    // Per-column objective-gradient flag (used in cb_read_matrix)
    const casadi_int* f_row = gradf_sp_.row();
    gradf_col_flag_.assign(nx_, false);
    gradf_col_to_nz_.assign(nx_, -1);
    for (casadi_int k = 0; k < gradf_sp_.nnz(); ++k) {
      gradf_col_flag_[f_row[k]] = true;
      gradf_col_to_nz_[f_row[k]] = k;
    }

    // Detect constant (linear) objective gradient entries using second-order sparsity
    {
      Sparsity dgradf_dx = gradf_fcn.sparsity_jac(0, 1, true);
      gradf_nlflag_.assign(gradf_sp_.nnz(), 0);
      const casadi_int* dg_row = dgradf_dx.row();
      for (casadi_int el = 0; el < dgradf_dx.nnz(); ++el) {
        gradf_nlflag_[dg_row[el]] = 1;
      }
      has_linear_gradf_ = std::any_of(gradf_nlflag_.begin(), gradf_nlflag_.end(),
                                       [](int f) { return f == 0; });
    }

    refine_nlflags_with_hessian();

    const casadi_int* g_colind = jacg_sp_.colind();
    const casadi_int* g_row = jacg_sp_.row();

    // Build row-indexed CSR structure over nonlinear entries for fast
    // Jacobian scatter in cb_fd_eval
    casadi_int nnz_g = jacg_sp_.nnz();
    jacg_rowstart_.assign(ng_ + 1, 0);
    for (casadi_int el = 0; el < nnz_g; ++el)
      if (jacg_nlflag_[el]) jacg_rowstart_[g_row[el] + 1]++;
    for (int r = 0; r < ng_; ++r)
      jacg_rowstart_[r + 1] += jacg_rowstart_[r];
    jacg_nzidx_.resize(jacg_rowstart_[ng_]);
    jacg_col_.resize(jacg_rowstart_[ng_]);
    std::vector<int> fill_pos(ng_, 0);
    for (int c = 0; c < nx_; ++c) {
      for (casadi_int el = g_colind[c]; el < g_colind[c+1]; ++el) {
        if (!jacg_nlflag_[el]) continue;
        int r = static_cast<int>(g_row[el]);
        casadi_assert(fill_pos[r] < jacg_rowstart_[r + 1] - jacg_rowstart_[r],
                      "CSR fill overflow for row r - count/fill pass mismatch in jacg_rowstart_");
        int pos = jacg_rowstart_[r] + fill_pos[r]++;
        jacg_nzidx_[pos] = static_cast<int>(el);
        jacg_col_[pos]   = c;
      }
    }

    build_tapes();
  }

  bool ConoptInterface::build_tape(const Function& f, ConoptTape& t) {
    t = ConoptTape();
    if (!f.is_a("SXFunction", false)) return false;
    const SXFunction* sxf = static_cast<const SXFunction*>(f.get());
    if (sxf->has_free()) return false;
    const std::vector<ScalarAtomic>& alg = sxf->algorithm_;
    casadi_int n = alg.size();
    if (n >= std::numeric_limits<int>::max()) return false;

    t.instr.resize(n);
    t.dep1.assign(n, -1);
    t.dep2.assign(n, -1);
    t.out_instr.resize(f.n_out());
    for (casadi_int i = 0; i < f.n_out(); ++i) t.out_instr[i].assign(f.nnz_out(i), -1);
    t.sz_w = sxf->worksize_;

    // Last instruction to write each work vector slot. The work vector reuses
    // slots (live variables), so an operand's producer is the most recent
    // writer of its slot; operands are resolved before the result is recorded
    // since an instruction may write the slot it reads.
    std::vector<int> writer(t.sz_w, -1);
    for (casadi_int k = 0; k < n; ++k) {
      const ScalarAtomic& e = alg[k];
      ConoptInstr& c = t.instr[k];
      c.op = e.op;
      c.i0 = e.i0;
      c.i1 = 0;
      c.i2 = 0;
      c.d = 0;
      switch (e.op) {
        case OP_CALL:
          return false;
        case OP_CONST:
          c.d = e.d;
          writer[e.i0] = static_cast<int>(k);
          break;
        case OP_INPUT:
          c.i1 = e.i1;
          c.i2 = e.i2;
          writer[e.i0] = static_cast<int>(k);
          break;
        case OP_OUTPUT:
          c.i1 = e.i1;
          c.i2 = e.i2;
          t.dep1[k] = writer[e.i1];
          t.out_instr[e.i0][e.i2] = static_cast<int>(k);
          break;
        default:
          // Unary operations have i2 == i1 (see SXFunction::init)
          c.i1 = e.i1;
          c.i2 = e.i2;
          t.dep1[k] = writer[e.i1];
          if (casadi_math<double>::ndeps(e.op) == 2) t.dep2[k] = writer[e.i2];
          writer[e.i0] = static_cast<int>(k);
      }
    }
    for (auto&& o : t.out_instr)
      if (std::any_of(o.begin(), o.end(), [](int k) { return k < 0; })) return false;
    return true;
  }

  void ConoptInterface::build_tapes() {
    has_tape_ = false;
    tape_fg_ = ConoptTape();
    tape_fg_jac_ = ConoptTape();
    if (!subset_eval_) return;
    // Row slots index g's nonzeros directly, as elsewhere in this interface
    has_tape_ = get_function("nlp_fg").nnz_out(1) == ng_ &&
                build_tape(get_function("nlp_fg"), tape_fg_) &&
                build_tape(get_function("nlp_fg_jac"), tape_fg_jac_);
    if (!has_tape_) {
      tape_fg_ = ConoptTape();
      tape_fg_jac_ = ConoptTape();
      if (verbose_ || debug_)
        casadi_message("CONOPT: row-subset evaluation disabled (evaluation functions are "
                       "not SX without calls); evaluating all rows.");
    }
  }

  void ConoptInterface::refine_nlflags_with_hessian() {
    // The second-order sparsity of jac:g:x / grad:f:x can be conservative, e.g.
    // through MX Function calls, where every input of the call gets flagged.
    // nlp_hess_l's structure is the union over all rows (symbolic multipliers),
    // so a column absent from it has zero second derivatives everywhere: all of
    // its Jacobian/gradient entries are constant.
    if (!exact_hessian_) return;
    std::vector<bool> in_hess(nx_, false);
    const casadi_int* h_colind = hesslag_sp_.colind();
    const casadi_int* h_row = hesslag_sp_.row();
    for (casadi_int c = 0; c < nx_; ++c) {
      for (casadi_int el = h_colind[c]; el < h_colind[c+1]; ++el) {
        in_hess[c] = true;
        in_hess[h_row[el]] = true;
      }
    }

    const casadi_int* g_colind = jacg_sp_.colind();
    for (casadi_int c = 0; c < nx_; ++c) {
      if (in_hess[c]) continue;
      for (casadi_int el = g_colind[c]; el < g_colind[c+1]; ++el) jacg_nlflag_[el] = 0;
    }
    const casadi_int* f_row = gradf_sp_.row();
    for (casadi_int k = 0; k < gradf_sp_.nnz(); ++k) {
      if (!in_hess[f_row[k]]) gradf_nlflag_[k] = 0;
    }

    has_linear_jac_ = std::any_of(jacg_nlflag_.begin(), jacg_nlflag_.end(),
                                   [](int f) { return f == 0; });
    has_linear_gradf_ = std::any_of(gradf_nlflag_.begin(), gradf_nlflag_.end(),
                                     [](int f) { return f == 0; });
  }

  // --- Serialization & Deserialization --- //
  ConoptInterface::ConoptInterface(DeserializingStream& s) : Nlpsol(s) {
    int version = s.version("ConoptInterface", 1, 2);
    s.unpack("ConoptInterface::exact_hessian", exact_hessian_);
    s.unpack("ConoptInterface::opts", opts_);
    s.unpack("ConoptInterface::gradf_sp", gradf_sp_);
    s.unpack("ConoptInterface::jacg_sp", jacg_sp_);
    s.unpack("ConoptInterface::hesslag_sp", hesslag_sp_);
    s.unpack("ConoptInterface::optfile", optfile_);
    s.unpack("ConoptInterface::warm_start", warm_start_);
    s.unpack("ConoptInterface::debug", debug_);
    subset_eval_ = true;
    if (version >= 2) s.unpack("ConoptInterface::subset_eval", subset_eval_);

    // Recompute linearity flags first (needed for CSR construction below)
    {
      Sparsity djac_dx = get_function("nlp_jac_g").sparsity_jac(0, 1, true);
      jacg_nlflag_.assign(jacg_sp_.nnz(), 0);
      const casadi_int* dj_row = djac_dx.row();
      for (casadi_int el = 0; el < djac_dx.nnz(); ++el)
        jacg_nlflag_[dj_row[el]] = 1;
      has_linear_jac_ = std::any_of(jacg_nlflag_.begin(), jacg_nlflag_.end(),
                                     [](int f) { return f == 0; });
    }

    // Rebuild derived arrays from the serialized sparsities
    const casadi_int* f_row = gradf_sp_.row();
    gradf_col_flag_.assign(nx_, false);
    gradf_col_to_nz_.assign(nx_, -1);
    for (casadi_int k = 0; k < gradf_sp_.nnz(); ++k) {
      gradf_col_flag_[f_row[k]] = true;
      gradf_col_to_nz_[f_row[k]] = k;
    }

    {
      Sparsity dgradf_dx = get_function("nlp_grad_f").sparsity_jac(0, 1, true);
      gradf_nlflag_.assign(gradf_sp_.nnz(), 0);
      const casadi_int* dg_row = dgradf_dx.row();
      for (casadi_int el = 0; el < dgradf_dx.nnz(); ++el)
        gradf_nlflag_[dg_row[el]] = 1;
      has_linear_gradf_ = std::any_of(gradf_nlflag_.begin(), gradf_nlflag_.end(),
                                       [](int f) { return f == 0; });
    }

    refine_nlflags_with_hessian();

    const casadi_int* g_colind = jacg_sp_.colind();
    const casadi_int* g_row = jacg_sp_.row();
    casadi_int nnz_g = jacg_sp_.nnz();
    jacg_rowstart_.assign(ng_ + 1, 0);
    for (casadi_int el = 0; el < nnz_g; ++el)
      if (jacg_nlflag_[el]) jacg_rowstart_[g_row[el] + 1]++;
    for (int r = 0; r < ng_; ++r)
      jacg_rowstart_[r + 1] += jacg_rowstart_[r];
    jacg_nzidx_.resize(jacg_rowstart_[ng_]);
    jacg_col_.resize(jacg_rowstart_[ng_]);
    std::vector<int> fill_pos(ng_, 0);
    for (int c = 0; c < nx_; ++c)
      for (casadi_int el = g_colind[c]; el < g_colind[c+1]; ++el) {
        if (!jacg_nlflag_[el]) continue;
        int r = static_cast<int>(g_row[el]);
        casadi_assert(fill_pos[r] < jacg_rowstart_[r + 1] - jacg_rowstart_[r],
                      "CSR fill overflow for row r - count/fill pass mismatch in jacg_rowstart_");
        int pos = jacg_rowstart_[r] + fill_pos[r]++;
        jacg_nzidx_[pos] = static_cast<int>(el);
        jacg_col_[pos]   = c;
      }

    build_tapes();
  }

  void ConoptInterface::serialize_body(SerializingStream &s) const {
    Nlpsol::serialize_body(s);
    s.version("ConoptInterface", 2);
    s.pack("ConoptInterface::exact_hessian", exact_hessian_);
    s.pack("ConoptInterface::opts", opts_);
    s.pack("ConoptInterface::gradf_sp", gradf_sp_);
    s.pack("ConoptInterface::jacg_sp", jacg_sp_);
    s.pack("ConoptInterface::hesslag_sp", hesslag_sp_);
    s.pack("ConoptInterface::optfile", optfile_);
    s.pack("ConoptInterface::warm_start", warm_start_);
    s.pack("ConoptInterface::debug", debug_);
    s.pack("ConoptInterface::subset_eval", subset_eval_);
    // Derived arrays (gradf_col_flag_, CSR, tapes) are rebuilt on deserialization
  }

  ConoptMemory::ConoptMemory(const ConoptInterface& interface)
      : self(interface), NlpsolMemory(), cntvect(nullptr),
        modsta(ConoptModelStatus::Unset), solsta(ConoptSolverStatus::Unset),
        iter(0), return_status("Unset"),
        cache_valid(false), cache_valid_jac(false), nan_encountered(false),
        ng_expanded(0), numnz_expanded(0) {}

  ConoptMemory::~ConoptMemory() {
    if (cntvect) COI_Free(&cntvect);
  }

  void ConoptMemory::bump_stamp() {
    if (++cur_stamp == 0) {
      // Wrapped around: clear so that no stale entry matches the new stamp
      std::fill(stamp_fun.begin(), stamp_fun.end(), 0);
      std::fill(stamp_jac.begin(), stamp_jac.end(), 0);
      cur_stamp = 1;
    }
  }

  void ConoptInterface::free_mem(void* mem) const { delete static_cast<ConoptMemory*>(mem); }

  int ConoptInterface::init_mem(void* mem) const {
    if (Nlpsol::init_mem(mem)) return 1;
    auto m = static_cast<ConoptMemory*>(mem);

    m->cached_x.resize(nx_, 0.0);
    m->cached_grad_f.resize(gradf_sp_.nnz(), 0.0);
    m->cached_g.resize(ng_, 0.0);
    m->cached_jac_g.resize(jacg_sp_.nnz(), 0.0);
    m->casadi_to_conopt_lb_row.resize(ng_);
    m->casadi_to_conopt_ub_row.assign(ng_, -1);
    m->hess_lam_g_.resize(ng_, 0.0);
    m->row_const_.assign(ng_, 0.0);
    m->row_nnz.assign(ng_, 0);
    // Every CasADi row contributes at least one entry (range constraints add a
    // second), so ng_ is a guaranteed lower bound on the final size — reserving
    // less would force a reallocation on essentially every solve. solve() grows
    // these vectors further on demand only for the range-constraint rows.
    casadi_int initial_row_reserve = ng_;
    m->conopt_to_casadi.reserve(initial_row_reserve);
    m->conopt_type.reserve(initial_row_reserve);
    m->conopt_rhs.reserve(initial_row_reserve);
    if (has_linear_jac_) {
      m->const_jac_vals.resize(jacg_sp_.nnz(), 0.0);
    }
    m->linear_at_x0.resize(ng_, 0.0);
    if (has_linear_gradf_)
      m->gradf_const_vals.resize(gradf_sp_.nnz(), 0.0);

    if (has_tape_) {
      m->stamp_fun.assign(ng_ + 1, 0);
      m->stamp_jac.assign(ng_ + 1, 0);
      m->queued.assign(ng_ + 1, 0);
      m->cur_stamp = 1;
      m->cur_queue = 0;
      m->pending_rows.reserve(ng_ + 1);
      m->tape_w.assign(std::max(tape_fg_.sz_w, tape_fg_jac_.sz_w), 0.0);
      m->tape_mark.assign(std::max(tape_fg_.instr.size(), tape_fg_jac_.instr.size()), 0);
    }

    if (COI_Create(&m->cntvect) != 0 || m->cntvect == nullptr) {
      casadi::uerr() << "CONOPT: COI_Create failed" << std::endl;
      return 1;
    }

    // License comes from the 'license' option if given, otherwise from
    // environment variables (see e.g. conopt/examples/setlicense.sh), rather
    // than being compiled in, so the same build works for any licensee.
    const char* lic_int1 = std::getenv("CONOPT_LICENSE_INT_1");
    const char* lic_int2 = std::getenv("CONOPT_LICENSE_INT_2");
    const char* lic_int3 = std::getenv("CONOPT_LICENSE_INT_3");
    const char* lic_text = std::getenv("CONOPT_LICENSE_TEXT");
    if (has_license_) {
      COIDEF_License(m->cntvect, license_int_[0], license_int_[1],
                      license_int_[2], license_text_.c_str());
    } else if (lic_int1 && lic_int2 && lic_int3 && lic_text) {
      COIDEF_License(m->cntvect, std::atoi(lic_int1), std::atoi(lic_int2),
                      std::atoi(lic_int3), lic_text);
    }

    if (warm_start_) COIDEF_IniStat(m->cntvect, 2);

    COIDEF_NumVar(m->cntvect, nx_);
    // NumCon, NumNz, NumNlNz are set in solve() because range-constraint expansion
    // can change them between calls.

    // CasADi allows variables that appear in no constraint and not in the
    // objective (or only with a zero constant gradient); these give empty columns.
    COIDEF_EmptyCol(m->cntvect, 1);

    COIDEF_ObjCon(m->cntvect, 0);
    COIDEF_OptDir(m->cntvect, -1);

    // Handle Options
    m->custom_options.clear();
    for (auto&& op : opts_) {
        // Explictly catch C API options defined in conopt.h
        if (op.first == "itlim") COIDEF_ItLim(m->cntvect, op.second.to_int());
        else if (op.first == "errlim") COIDEF_ErrLim(m->cntvect, op.second.to_int());
        else if (op.first == "reslim" || op.first == "timelim")
            COIDEF_ResLim(m->cntvect, op.second.to_double());
        else if (op.first == "maxheap") COIDEF_MaxHeap(m->cntvect, op.second.to_double());
        else if (op.second.is_string()) {
            casadi_warning("CONOPT option '" + op.first + "' is a string; string options cannot be "
                           "passed via the CONOPT option callback (no SVAL parameter). "
                           "Use the 'optfile' option instead.");
        } else {
            m->custom_options.push_back(op);
        }
    }
    if (!optfile_.empty()) COIDEF_Optfile(m->cntvect, optfile_.c_str());
    COIDEF_Option(m->cntvect, &ConoptInterface::cb_option);
    COIDEF_Progress(m->cntvect, &ConoptInterface::cb_progress);

    // Register Callbacks
    COIDEF_ReadMatrix(m->cntvect, &ConoptInterface::cb_read_matrix);
    COIDEF_FDEvalIni(m->cntvect, &ConoptInterface::cb_fdevalini);
    COIDEF_FDEval(m->cntvect, &ConoptInterface::cb_fd_eval);
    COIDEF_FDEvalEnd(m->cntvect, &ConoptInterface::cb_fdevalend);

    if (exact_hessian_ && hesslag_sp_.nnz() > 0) {
        COIDEF_NumHess(m->cntvect, hesslag_sp_.nnz());
        COIDEF_2DLagrStr(m->cntvect, &ConoptInterface::cb_2dlagrstr);
        COIDEF_2DLagrVal(m->cntvect, &ConoptInterface::cb_2dlagrval);
    }

    COIDEF_FVincLin(m->cntvect, 1);

    COIDEF_Status(m->cntvect, &ConoptInterface::cb_status);
    COIDEF_Solution(m->cntvect, &ConoptInterface::cb_solution);
    COIDEF_Message(m->cntvect, &ConoptInterface::cb_message);
    COIDEF_ErrMsg(m->cntvect, &ConoptInterface::cb_errmsg);
    COIDEF_UsrMem(m->cntvect, m);

    return 0;
  }

  void ConoptInterface::set_work(void* mem, const double**& arg, double**& res,
                                  casadi_int*& iw, double*& w) const {
    Nlpsol::set_work(mem, arg, res, iw, w);
  }

  // conopt_to_casadi/conopt_type/conopt_rhs always grow in lockstep, so a single
  // capacity check (on conopt_to_casadi) is enough to decide whether to grow all
  // three. Growth is by 0.25 of the CasADi rows not yet processed, rather than
  // jumping straight to the worst-case (all-range) size.
  void ConoptInterface::ensure_row_capacity(ConoptMemory* m, casadi_int remaining_rows) const {
    if (m->conopt_to_casadi.size() == m->conopt_to_casadi.capacity()) {
      casadi_int growth = std::max<casadi_int>(
          std::min<casadi_int>(remaining_rows, 10), remaining_rows / 4);
      casadi_int new_cap = m->conopt_to_casadi.capacity() + growth;
      m->conopt_to_casadi.reserve(new_cap);
      m->conopt_type.reserve(new_cap);
      m->conopt_rhs.reserve(new_cap);
    }
  }

  int ConoptInterface::solve(void* mem) const {
    auto m = static_cast<ConoptMemory*>(mem);
    m->cache_valid     = false;
    m->cache_valid_jac = false;
    // Parameters/bounds may differ from the previous solve, so a cache for the
    // same x is not reusable across solves.
    m->stored_fun      = false;
    m->stored_jac      = false;
    if (has_tape_) m->bump_stamp();
    m->n_eval_subset = m->n_eval_full = m->n_eval_reused = m->n_eval_on_demand = 0;
    m->subset_instr_frac = 0;
    m->cached_f        = 0.0;
    m->nan_encountered = false;
    m->modsta = ConoptModelStatus::Unset;
    m->solsta = ConoptSolverStatus::Unset;
    m->iter = 0;
    m->return_status = "Unset";

    // Build the per-solve constraint expansion (splits range constraints into two rows)
    std::fill(m->row_const_.begin(), m->row_const_.end(), 0.0);
    m->obj_const_lin_ = 0.0;
    m->conopt_to_casadi.clear();
    m->casadi_to_conopt_ub_row.assign(ng_, -1);
    m->conopt_type.clear();
    m->conopt_rhs.clear();

    // Compute total nnz per CasADi row (needed for range-constraint NZ duplication)
    std::fill(m->row_nnz.begin(), m->row_nnz.end(), 0);
    {
      const casadi_int* g_row_s = jacg_sp_.row();
      for (casadi_int el = 0; el < (casadi_int)jacg_sp_.nnz(); ++el)
        m->row_nnz[g_row_s[el]]++;
    }

    casadi_int ng_expanded = 0;
    // Objective gradient NZs are added to numnz after nlp_grad_f is evaluated
    // (so that gradf_const_vals is populated before we check for non-zero linear entries).
    casadi_int numnz = (casadi_int)jacg_sp_.nnz();

    for (casadi_int i = 0; i < ng_; ++i) {
      double lbg = m->d_nlp.lbz[nx_ + i];
      double ubg = m->d_nlp.ubz[nx_ + i];
      bool is_range = !std::isinf(lbg) && !std::isinf(ubg) && lbg != ubg;

      // CONOPT row 0 is reserved for the objective, so constraint rows start at 1
      // (arrays are still plain 0-based C arrays; only the row *numbering* is offset).
      m->casadi_to_conopt_lb_row[i] = static_cast<int>(ng_expanded + 1);
      ensure_row_capacity(m, ng_ - i);
      m->conopt_to_casadi.push_back(static_cast<int>(i));
      if (lbg == ubg) {
        m->conopt_type.push_back(ConoptRowType::Equality);  m->conopt_rhs.push_back(lbg);
      } else if (!std::isinf(lbg) && std::isinf(ubg)) {
        m->conopt_type.push_back(ConoptRowType::GreaterEqual);  m->conopt_rhs.push_back(lbg);
      } else if (std::isinf(lbg) && !std::isinf(ubg)) {
        m->conopt_type.push_back(ConoptRowType::LessEqual);  m->conopt_rhs.push_back(ubg);
      } else if (std::isinf(lbg) && std::isinf(ubg)) {
        m->conopt_type.push_back(ConoptRowType::Free);  m->conopt_rhs.push_back(0.0);
      } else {
        // range: >= row
        m->conopt_type.push_back(ConoptRowType::GreaterEqual);
        m->conopt_rhs.push_back(lbg);
      }
      ng_expanded++;

      if (is_range) {
        m->casadi_to_conopt_ub_row[i] = static_cast<int>(ng_expanded + 1);
        ensure_row_capacity(m, ng_ - i);
        m->conopt_to_casadi.push_back(static_cast<int>(i)); // <= row
        m->conopt_type.push_back(ConoptRowType::LessEqual);
        m->conopt_rhs.push_back(ubg);
        ng_expanded++;
        numnz += m->row_nnz[i];
      }
    }
    casadi_assert(ng_expanded <= std::numeric_limits<int>::max(), "ng_expanded overflows int");
    m->ng_expanded = static_cast<int>(ng_expanded);

    // Empty Jacobian rows also need their constant terms moved into the RHS.
    bool has_affine_g = has_linear_jac_ ||
        std::any_of(m->row_nnz.begin(), m->row_nnz.end(), [](int n) { return n == 0; });
    if (has_affine_g) {
      m->arg[0] = m->d_nlp.z;
      m->arg[1] = m->d_nlp.p;
      m->res[0] = m->cached_g.data();
      m->res[1] = has_linear_jac_ ? m->const_jac_vals.data() : nullptr;
      try {
        if (calc_function(m, "nlp_jac_g")) {
          m->success = false;
          m->unified_return_status = SOLVER_RET_NAN;
          m->return_status = "Initial evaluation failed";
          return 0;
        }
      } catch (std::exception& ex) {
        casadi::uerr() << "CONOPT: initial evaluation failed: " << ex.what() << std::endl;
        return 1;
      } catch (...) {
        casadi::uerr() << "CONOPT: initial evaluation failed (unknown exception)" << std::endl;
        return 1;
      }
    }

    // CONOPT evaluates affine rows internally; absorb their constants into the RHS.
    if (has_affine_g) {
      const casadi_int* g_colind_c = jacg_sp_.colind();
      const casadi_int* g_row_c    = jacg_sp_.row();

      // Accumulate the linear part of G at x0 per row: sum_j a_j * x0_j
      std::fill(m->linear_at_x0.begin(), m->linear_at_x0.end(), 0.0);
      for (int c = 0; c < nx_; ++c) {
        for (casadi_int el = g_colind_c[c]; el < g_colind_c[c + 1]; ++el) {
          if (jacg_nlflag_[el] == 0)
            m->linear_at_x0[g_row_c[el]] += m->const_jac_vals[el] * m->d_nlp.z[c];
        }
      }

      for (int ci = 0; ci < ng_; ++ci) {
        // Only adjust fully linear rows (no nonlinear Jacobian entries)
        if (jacg_rowstart_[ci + 1] != jacg_rowstart_[ci]) continue;
        double constant = m->cached_g[ci] - m->linear_at_x0[ci];
        if (std::abs(constant) < 1e-14) continue;
        m->row_const_[ci] = constant;
        int lb_row = m->casadi_to_conopt_lb_row[ci];
        m->conopt_rhs[lb_row - 1] -= constant;
        int ub_row = m->casadi_to_conopt_ub_row[ci];
        if (ub_row >= 0) m->conopt_rhs[ub_row - 1] -= constant;
      }
    }

    // Evaluate objective gradient at initial point for constant (linear) entries,
    // or to capture the function value when the gradient is structurally empty.
    m->obj_const_ = std::numeric_limits<double>::quiet_NaN();
    if (has_linear_gradf_ || gradf_sp_.nnz() == 0) {
      m->arg[0] = m->d_nlp.z;
      m->arg[1] = m->d_nlp.p;
      m->res[0] = &m->cached_f;
      m->res[1] = has_linear_gradf_ ? m->gradf_const_vals.data() : nullptr;
      try {
        if (calc_function(m, "nlp_grad_f")) {
          m->success = false;
          m->unified_return_status = SOLVER_RET_NAN;
          m->return_status = "Initial evaluation failed";
          return 0;
        }
      } catch (std::exception& ex) {
        casadi::uerr() << "CONOPT: initial evaluation failed: " << ex.what() << std::endl;
        return 1;
      } catch (...) {
        casadi::uerr() << "CONOPT: initial evaluation failed (unknown exception)" << std::endl;
        return 1;
      }
    }

    // Detect constant objective at solve time: no nonlinear gradient entries and
    // all linear-gradient values are zero (objective has no x-dependence).
    // Switch to feasibility mode so CONOPT doesn't report OBJVAL=0 for an empty row.
    {
      bool has_nl_gradf = std::any_of(gradf_nlflag_.begin(), gradf_nlflag_.end(),
                                       [](int f) { return f == 1; });
      bool all_const_zero = has_linear_gradf_ &&
          std::all_of(m->gradf_const_vals.begin(), m->gradf_const_vals.end(),
                      [](double v) { return v == 0.0; });
      if (!has_nl_gradf && (gradf_sp_.nnz() == 0 || all_const_zero)) {
        m->obj_const_ = m->cached_f;
        COIDEF_OptDir(m->cntvect, 0);
      } else {
        if (!has_nl_gradf) {
          // CONOPT omits the affine objective's constant term; recover it at x0.
          const casadi_int* f_row_c = gradf_sp_.row();
          double lin_at_x0 = 0.0;
          for (casadi_int k = 0; k < gradf_sp_.nnz(); ++k)
            lin_at_x0 += m->gradf_const_vals[k] * m->d_nlp.z[f_row_c[k]];
          m->obj_const_lin_ = m->cached_f - lin_at_x0;
        }
        // cntvect persists across solves, so explicitly reset OptDir in case a
        // prior solve on this instance (e.g. with different parameters) hit the
        // constant-objective branch above and left it at 0.
        COIDEF_OptDir(m->cntvect, -1);
      }
    }

    // Count nonlinear NZ: nonlinear objective gradient entries + nonlinear constraint entries
    const casadi_int* g_colind_s = jacg_sp_.colind();
    const casadi_int* g_row_s    = jacg_sp_.row();
    casadi_int num_nl_nz = 0;
    for (casadi_int k = 0; k < (casadi_int)gradf_sp_.nnz(); ++k)
      if (gradf_nlflag_[k]) num_nl_nz++;
    for (casadi_int c = 0; c < nx_; ++c) {
      for (casadi_int el = g_colind_s[c]; el < g_colind_s[c+1]; ++el) {
        if (jacg_nlflag_[el]) {
          int ci = static_cast<int>(g_row_s[el]);
          num_nl_nz += (m->casadi_to_conopt_ub_row[ci] >= 0) ? 2 : 1;
        }
      }
    }

    casadi_assert(num_nl_nz <= std::numeric_limits<int>::max(), "num_nl_nz overflows int");

    // Count objective gradient NZs now that gradf_const_vals has been populated.
    // Nonlinear entries are always included; linear (constant) entries only when
    // non-zero — a zero constant gradient contributes nothing to the objective row
    // and must not occupy a slot in the CONOPT matrix structure.
    {
      casadi_int numnz_f = 0;
      for (casadi_int k = 0; k < (casadi_int)gradf_sp_.nnz(); ++k) {
        if (gradf_nlflag_[k] == 1 ||
            (has_linear_gradf_ && m->gradf_const_vals[k] != 0.0))
          numnz_f++;
      }
      numnz += numnz_f;
    }
    casadi_assert(numnz <= std::numeric_limits<int>::max(), "numnz overflows int");
    m->numnz_expanded = static_cast<int>(numnz);

    COIDEF_NumCon(m->cntvect, static_cast<int>(ng_expanded + 1));
    COIDEF_NumNz(m->cntvect, static_cast<int>(numnz));
    COIDEF_NumNlNz(m->cntvect, static_cast<int>(num_nl_nz));

    int ret = COI_Solve(m->cntvect);

    // Restore the affine objective's constant term (cb_status only saw grad'x).
    m->d_nlp.objective += m->obj_const_lin_;

    // Restore constant objective value when CONOPT ran in feasibility mode.
    if (!std::isnan(m->obj_const_)) m->d_nlp.objective = m->obj_const_;

    if (ret != 0) {
      m->success = false;
      m->unified_return_status = m->nan_encountered ? SOLVER_RET_NAN : SOLVER_RET_UNKNOWN;
      return 0;
    }

    m->success = !m->nan_encountered &&
                 (m->modsta == ConoptModelStatus::Optimal ||
                  m->modsta == ConoptModelStatus::LocallyOptimal) &&
                 m->solsta == ConoptSolverStatus::NormalCompletion;

    if (m->nan_encountered || m->solsta == ConoptSolverStatus::EvalErrorLimit) {
      m->unified_return_status = SOLVER_RET_NAN;
    } else if (m->success) {
      m->unified_return_status = SOLVER_RET_SUCCESS;
    } else {
      if (m->solsta == ConoptSolverStatus::IterationLimit ||
          m->solsta == ConoptSolverStatus::TimeLimit ||
          m->solsta == ConoptSolverStatus::UserInterrupt ||
          m->solsta == ConoptSolverStatus::QuickModeTermination) {
        m->unified_return_status = SOLVER_RET_LIMITED;
      } else if (m->modsta == ConoptModelStatus::Infeasible ||
                 m->modsta == ConoptModelStatus::LocallyInfeasible) {
        m->unified_return_status = SOLVER_RET_INFEASIBLE;
      }
    }
    return 0;
  }

  Dict ConoptInterface::get_stats(void* mem) const {
    Dict stats = Nlpsol::get_stats(mem);
    auto m = static_cast<ConoptMemory*>(mem);
    stats["return_status"] = m->return_status;
    stats["modsta"] = static_cast<int>(m->modsta);
    stats["solsta"] = static_cast<int>(m->solsta);
    stats["iter_count"] = m->iter;
    stats["subset_eval"] = has_tape_;
    stats["n_eval_subset"] = m->n_eval_subset;
    stats["n_eval_full"] = m->n_eval_full;
    stats["n_eval_reused"] = m->n_eval_reused;
    stats["n_eval_on_demand"] = m->n_eval_on_demand;
    stats["subset_instr_frac"] = m->n_eval_subset > 0 ?
        m->subset_instr_frac / static_cast<double>(m->n_eval_subset) : 0.0;
    return stats;
  }

  // --- Dynamic Option Callback --- //
  int COI_CALLCONV ConoptInterface::cb_option(int NCALL, double* RVAL, int* IVAL,
                                               int* LVAL, char* NAME, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);

    if (NCALL >= static_cast<int>(m->custom_options.size())) {
        NAME[0] = '\0';
        return 0;
    }

    auto& opt = m->custom_options[NCALL];
    // NAME has a fixed-size buffer; init() validates the option-name length.
    size_t name_len = std::min(opt.first.size(), conopt_max_option_name);
    std::memcpy(NAME, opt.first.c_str(), name_len);
    NAME[name_len] = '\0';

    if (opt.second.is_double()) {
        *RVAL = opt.second.to_double();
        if (m->self.debug_)
            casadi::uout() << "CONOPT option: " << opt.first << " = " << *RVAL << std::endl;
    } else if (opt.second.is_int()) {
        *IVAL = opt.second.to_int();
        if (m->self.debug_)
            casadi::uout() << "CONOPT option: " << opt.first << " = " << *IVAL << std::endl;
    } else if (opt.second.is_bool()) {
        *LVAL = opt.second.to_bool() ? 1 : 0;
        if (m->self.debug_) {
            casadi::uout() << "CONOPT option: " << opt.first << " = "
                           << (opt.second.to_bool() ? "true" : "false") << std::endl;
        }
    } else if (opt.second.is_string()) {
        // init_mem filters out all string options before they enter custom_options
        // (with a casadi_warning directing the user to 'optfile'). Reaching this
        // branch means custom_options was populated externally in a way that
        // bypasses that filter, which is a programming error. Setting NAME[0]='\0'
        // here would terminate the entire option enumeration, silently dropping
        // all subsequent entries — so we assert rather than pretend to skip.
        casadi_error("CONOPT option '" + opt.first + "' is a string type in cb_option. "
                     "String options cannot be passed via the CONOPT option callback "
                     "(COI_OPTION_t has no SVAL parameter). Use the 'optfile' option "
                     "instead. The init_mem filter should have removed this option "
                     "before it reached custom_options; reaching this branch is a "
                     "programming error.");
    } else {
        // Similarly, an option of unknown type must never reach this point.
        casadi_error("CONOPT option '" + opt.first + "' has an unknown GenericType in "
                     "cb_option. Only double, int, and bool options are valid here. "
                     "Reaching this branch is a programming error.");
    }

    return 0;
  }

  // --- Progress / Interrupt Callback --- //
  int COI_CALLCONV ConoptInterface::cb_progress(int LEN_INT, const int INTX[], int LEN_RL,
                                                 const double RL[], const double X[],
                                                 void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    const ConoptInterface& self = m->self;

    if (!self.fcallback_.is_null()) {
        int phase = (LEN_INT > 1) ? INTX[1] : -1;

        double obj_val = m->cached_f;  // best available approximation for early phases
        if (LEN_RL > 1 && phase >= 3) {
            obj_val = RL[1];            // CONOPT-reported value once available
        }

        std::fill_n(m->arg, self.fcallback_.n_in(), nullptr);
        m->arg[NLPSOL_X] = X;
        m->arg[NLPSOL_F] = &obj_val;

        std::fill_n(m->res, self.fcallback_.n_out(), nullptr);
        double ret_double = 0;
        m->res[0] = &ret_double;

        try {
            self.fcallback_(m->arg, m->res, m->iw, m->w, 0);
            if (ret_double != 0.0) return 1;
        } catch (KeyboardInterruptException& ex) {
            return 1;
        } catch (std::exception& ex) {
            casadi_warning(std::string("intermediate_callback: ") + ex.what());
            if (!self.iteration_callback_ignore_errors_) return 1;
        }
    }
    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_read_matrix(double LOWER[], double CURR[],
                                                    double UPPER[], int VSTA[], int TYPEX[],
                                                    double RHS[], int ESTA[], int COLSTA[],
                                                    int ROWNO[], double VALUE[],
                                                    int NLFLAG[], int NUMVAR, int NUMCON,
                                                    int NUMNZ, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    const ConoptInterface& self = m->self;

    // CONOPT is expected to echo back exactly what we told it via COIDEF_NumVar/
    // COIDEF_NumCon in init_mem()/solve() — a mismatch means the interface's own
    // bookkeeping (nx_, ng_expanded) is desynced from what CONOPT is using.
    casadi_assert(NUMVAR == self.nx_, "cb_read_matrix: NUMVAR != nx_");
    casadi_assert(NUMCON == m->ng_expanded + 1, "cb_read_matrix: NUMCON != ng_expanded + 1");

    // Variable bounds and initial point (clamped to [lb, ub])
    for (int i = 0; i < NUMVAR; ++i) {
      double lb = m->d_nlp.lbz[i];
      double ub = m->d_nlp.ubz[i];
      if (!std::isinf(lb)) LOWER[i] = lb;
      if (!std::isinf(ub)) UPPER[i] = ub;
      double x0 = m->d_nlp.z[i];
      if (!std::isinf(ub)) x0 = std::min(x0, ub);
      if (!std::isinf(lb)) x0 = std::max(x0, lb);
      CURR[i] = x0;
    }

    // Constraint types and RHS (row 0 = objective, rows 1..ng_expanded = constraints).
    // conopt_type/conopt_rhs come from solve()'s range-constraint expansion —
    // they depend on this call's numeric lbg/ubg, not just problem structure.
    TYPEX[0] = static_cast<int>(ConoptRowType::Free);
    RHS[0]   = 0.0;
    for (int r = 0; r < m->ng_expanded; ++r) {
      TYPEX[r + 1] = static_cast<int>(m->conopt_type[r]);
      RHS[r + 1]   = m->conopt_rhs[r];
    }

    if (self.warm_start_) {
      // conopt_to_casadi/casadi_to_conopt_lb_row/ub_row (from solve()) map each
      // expanded CONOPT row back to its CasADi constraint, to look up lam there.
      ESTA[0] = static_cast<int>(ConoptBasisStatus::SuperBasic);  // objective row always superbasic

      // lam is the prior solve's dual solution, length nx_+ng_. lam[0..nx_-1] are
      // variable-bound multipliers, lam[nx_..] are constraint multipliers.
      const double* lam = m->d_nlp.lam;
      bool all_zero = std::all_of(lam, lam + self.nx_ + self.ng_,
                                  [](double v) { return v == 0.0; });
      if (!all_zero) {
        for (int i = 0; i < NUMVAR; ++i) {
          double lbi = m->d_nlp.lbz[i];
          double ubi = m->d_nlp.ubz[i];
          double xi  = m->d_nlp.z[i];
          double li  = lam[i];  // nonzero multiplier => that bound is active, so nonbasic
          if (!std::isinf(lbi) && std::fabs(xi - lbi) < 1e-8 && li <= 0.0)
            VSTA[i] = static_cast<int>(ConoptBasisStatus::AtLower);
          else if (!std::isinf(ubi) && std::fabs(xi - ubi) < 1e-8 && li >= 0.0)
            VSTA[i] = static_cast<int>(ConoptBasisStatus::AtUpper);
          else
            VSTA[i] = static_cast<int>(ConoptBasisStatus::Basic);
        }

        for (int r = 0; r < m->ng_expanded; ++r) {
          int ci        = m->conopt_to_casadi[r];
          // this constraint's multiplier, same sign logic as above
          double lam_ci = lam[NUMVAR + ci];
          int row1      = m->casadi_to_conopt_lb_row[ci];
          int row2      = m->casadi_to_conopt_ub_row[ci];
          if (r + 1 == row1) {
            ConoptRowType ctype = m->conopt_type[r];  // type of this CONOPT expanded row
            if (ctype == ConoptRowType::Equality) {       // equality: both sides, just mark basic
              ESTA[r + 1] = static_cast<int>(ConoptBasisStatus::Basic);
            } else if (ctype == ConoptRowType::GreaterEqual) {  // >= row: active when lam_ci < 0
              ESTA[r + 1] = static_cast<int>(
                  (lam_ci < 0.0) ? ConoptBasisStatus::AtLower : ConoptBasisStatus::Basic);
            } else if (ctype == ConoptRowType::LessEqual) {
              // <= row (pure <= stored as lb_row): active when lam_ci > 0
              ESTA[r + 1] = static_cast<int>(
                  (lam_ci > 0.0) ? ConoptBasisStatus::AtUpper : ConoptBasisStatus::Basic);
            } else {                        // free row: superbasic
              ESTA[r + 1] = static_cast<int>(ConoptBasisStatus::SuperBasic);
            }
          } else if (r + 1 == row2) {
            // range ub side is active when lam_ci > 0
            ESTA[r + 1] = static_cast<int>(
                (lam_ci > 0.0) ? ConoptBasisStatus::AtUpper : ConoptBasisStatus::Basic);
          } else {
            ESTA[r + 1] = static_cast<int>(ConoptBasisStatus::Basic);
          }
        }
      }
    }

    // Jacobian structure — built live from jacg_sp_ with range-row duplication.
    // Column/sparsity/linearity data (jacg_sp_, gradf_col_flag_, gradf_col_to_nz_,
    // gradf_nlflag_, jacg_nlflag_) is structural, fixed since init(). Row mapping
    // (casadi_to_conopt_lb_row/ub_row) and numeric values (gradf_const_vals,
    // const_jac_vals) are this solve()'s numeric data.
    const casadi_int* g_colind = self.jacg_sp_.colind();
    const casadi_int* g_row    = self.jacg_sp_.row();
    int nz = 0;
    for (int c = 0; c < NUMVAR; ++c) {
      COLSTA[c] = nz;
      if (self.gradf_col_flag_[c]) {
        // Objective-gradient entry in column c, if any (row 0 = objective).
        casadi_int k = self.gradf_col_to_nz_[c];
        if (self.gradf_nlflag_[k] == 1) {
          ROWNO[nz]  = 0;
          NLFLAG[nz] = 1;
          nz++;
        } else if (m->gradf_const_vals[k] != 0.0) {
          ROWNO[nz]  = 0;
          NLFLAG[nz] = 0;
          VALUE[nz]  = m->gradf_const_vals[k];
          nz++;
        }
      }
      // Constraint-Jacobian entries in column c; duplicated onto the ub row too
      // when this CasADi row was split into a range constraint by solve().
      for (casadi_int el = g_colind[c]; el < g_colind[c+1]; ++el) {
        int ci = static_cast<int>(g_row[el]);
        int nlflag = self.jacg_nlflag_[el];
        ROWNO[nz]  = m->casadi_to_conopt_lb_row[ci];
        NLFLAG[nz] = nlflag;
        if (nlflag == 0) VALUE[nz] = m->const_jac_vals[el];
        nz++;
        if (m->casadi_to_conopt_ub_row[ci] >= 0) {
          ROWNO[nz]  = m->casadi_to_conopt_ub_row[ci];
          NLFLAG[nz] = nlflag;
          if (nlflag == 0) VALUE[nz] = m->const_jac_vals[el];
          nz++;
        }
      }
    }
    COLSTA[NUMVAR] = nz;

    casadi_assert(nz == NUMNZ,
                   "cb_read_matrix: nz != NUMNZ - Jacobian nonzero count mismatch "
                   "between solve()'s numnz computation and this callback's writes");

    return 0;
  }

  // Marks the row slots in m->pending_rows (or all slots) as cached at cached_x.
  static void mark_rows_cached(ConoptMemory* m, bool need_jac, bool all) {
    if (all) {
      std::fill(m->stamp_fun.begin(), m->stamp_fun.end(), m->cur_stamp);
      if (need_jac) std::fill(m->stamp_jac.begin(), m->stamp_jac.end(), m->cur_stamp);
    } else {
      for (int r : m->pending_rows) {
        m->stamp_fun[r] = m->cur_stamp;
        if (need_jac) m->stamp_jac[r] = m->cur_stamp;
      }
    }
  }

  int ConoptInterface::eval_full(ConoptMemory* m, bool need_jac) const {
    m->arg[0] = m->cached_x.data();
    m->arg[1] = m->d_nlp.p;
    m->res[0] = &m->cached_f;
    if (need_jac) {
      m->res[1] = m->cached_grad_f.data();
      m->res[2] = m->cached_g.data();
      m->res[3] = m->cached_jac_g.data();
      return calc_function(m, "nlp_fg_jac");
    } else {
      m->res[1] = m->cached_g.data();
      return calc_function(m, "nlp_fg");
    }
  }

  int ConoptInterface::eval_rows(ConoptMemory* m, bool need_jac) const {
    const ConoptTape& t = need_jac ? tape_fg_jac_ : tape_fg_;
    const size_t n_instr = t.instr.size();
    // Beyond this, the indirect loop outweighs what is skipped
    const size_t max_instr = n_instr / 2;

    // Row sets repeat, so reuse what was found for this set before
    std::sort(m->pending_rows.begin(), m->pending_rows.end());
    auto& sets = m->row_sets[need_jac ? 1 : 0];
    auto it = sets.find(m->pending_rows);
    if (it != sets.end() && it->second.heavy) return -1;

    // Timed and counted as the full evaluation function, for comparable stats
    ScopedTiming tic(m->thread_local_mem.at(0)->fstats.at(need_jac ? "nlp_fg_jac" : "nlp_fg"));
    InterruptHandler::check();

    const std::vector<int>* list = nullptr;
    if (it != sets.end()) {
      list = &it->second.instr;
    } else {
      // Output index of g in the tape's function
      const casadi_int g_out = need_jac ? 2 : 1;
      char* mark = m->tape_mark.data();
      std::vector<int>& stack = m->tape_stack;
      std::vector<int>& found = m->tape_list;
      stack.clear();
      found.clear();
      auto push = [&](int k) {
        if (k >= 0 && !mark[k]) {
          mark[k] = 1;
          stack.push_back(k);
        }
      };

      // Seed with the output instructions of the requested values and nonlinear
      // derivatives (linear entries are constants held by CONOPT).
      for (int r : m->pending_rows) {
        if (r == ng_) {
          for (int k : t.out_instr[0]) push(k);
          if (need_jac) {
            for (casadi_int k = 0; k < gradf_sp_.nnz(); ++k)
              if (gradf_nlflag_[k]) push(t.out_instr[1][k]);
          }
        } else {
          push(t.out_instr[g_out][r]);
          if (need_jac) {
            for (int k = jacg_rowstart_[r]; k < jacg_rowstart_[r + 1]; ++k)
              push(t.out_instr[3][jacg_nzidx_[k]]);
          }
        }
      }

      // Collect every instruction the outputs depend on; the cost is proportional
      // to the number collected, and stops once past max_instr.
      bool heavy = false;
      while (!stack.empty()) {
        if (found.size() > max_instr) {
          heavy = true;
          break;
        }
        int k = stack.back();
        stack.pop_back();
        found.push_back(k);
        push(t.dep1[k]);
        push(t.dep2[k]);
      }
      for (int k : found) mark[k] = 0;
      for (int k : stack) mark[k] = 0;

      // The tape is in topological order
      if (!heavy) std::sort(found.begin(), found.end());

      // Remember the outcome, bounding the memory held by stored lists
      size_t cap = std::max<size_t>(size_t(1) << 20, 4 * n_instr);
      if (m->row_set_ints + found.size() > cap) {
        m->row_sets[0].clear();
        m->row_sets[1].clear();
        m->row_set_ints = 0;
      }
      ConoptRowSet& rs = sets[m->pending_rows];
      rs.heavy = heavy;
      if (heavy) return -1;
      rs.instr = found;
      m->row_set_ints += found.size();
      list = &rs.instr;
    }

    m->n_eval_subset++;
    m->subset_instr_frac += static_cast<double>(list->size()) / static_cast<double>(n_instr);

    double* w = m->tape_w.data();
    const double* arg[2] = {m->cached_x.data(), m->d_nlp.p};
    double* res[4];
    res[0] = &m->cached_f;
    if (need_jac) {
      res[1] = m->cached_grad_f.data();
      res[2] = m->cached_g.data();
      res[3] = m->cached_jac_g.data();
    } else {
      res[1] = m->cached_g.data();
    }
    const ConoptInstr* instr = t.instr.data();
    for (int k : *list) {
      const ConoptInstr& e = instr[k];
      switch (e.op) {
        CASADI_MATH_FUN_BUILTIN(w[e.i1], w[e.i2], w[e.i0])

        case OP_CONST: w[e.i0] = e.d; break;
        case OP_INPUT: w[e.i0] = arg[e.i1] == nullptr ? 0 : arg[e.i1][e.i2]; break;
        case OP_OUTPUT: res[e.i0][e.i2] = w[e.i1]; break;
        default:
          casadi_error("Unknown operation " + str(e.op));
      }
    }

    // Make sure the computed entries are not NaN or Inf
    auto bad = [](double v) { return !std::isfinite(v); };
    for (int r : m->pending_rows) {
      bool fail = false;
      if (r == ng_) {
        fail = !t.out_instr[0].empty() && bad(m->cached_f);
        if (need_jac) {
          for (casadi_int k = 0; k < gradf_sp_.nnz() && !fail; ++k)
            fail = gradf_nlflag_[k] && bad(m->cached_grad_f[k]);
        }
      } else {
        fail = bad(m->cached_g[r]);
        if (need_jac) {
          for (int k = jacg_rowstart_[r]; k < jacg_rowstart_[r + 1] && !fail; ++k)
            fail = bad(m->cached_jac_g[jacg_nzidx_[k]]);
        }
      }
      if (fail) {
        if (debug_) {
          casadi::uout() << "Row-subset evaluation: NaN/Inf in "
                         << (r == ng_ ? std::string("objective") : "g[" + str(r) + "]")
                         << "\n";
        }
        return 1;
      }
    }
    return 0;
  }

  int ConoptInterface::fdevalini_subset(ConoptMemory* m, const double X[],
                                        const int ROWLIST[], int MODE, int LISTSIZE) const {
    const bool need_jac = (MODE != 1);

    if (std::memcmp(m->cached_x.data(), X, nx_ * sizeof(double)) != 0) {
      std::memcpy(m->cached_x.data(), X, nx_ * sizeof(double));
      m->bump_stamp();
    }

    // Row slots in ROWLIST not yet cached at this point, without duplicates
    // (both CONOPT rows of a range constraint map to the same slot).
    const std::vector<unsigned>& stamp = need_jac ? m->stamp_jac : m->stamp_fun;
    if (++m->cur_queue == 0) {
      std::fill(m->queued.begin(), m->queued.end(), 0);
      m->cur_queue = 1;
    }
    m->pending_rows.clear();
    bool all = false;
    for (int i = 0; i < LISTSIZE; ++i) {
      int row = ROWLIST[i];
      if (row < 0 || row > m->ng_expanded) {
        all = true;  // unexpected numbering: be safe
        break;
      }
      int r = row == 0 ? static_cast<int>(ng_) : m->conopt_to_casadi[row - 1];
      if (stamp[r] == m->cur_stamp || m->queued[r] == m->cur_queue) continue;
      m->queued[r] = m->cur_queue;
      m->pending_rows.push_back(r);
    }

    if (!all && m->pending_rows.empty()) {
      m->n_eval_reused++;
      if (debug_) casadi::uout() << "FDEvalIni: requested rows already evaluated at this point\n";
      return 0;
    }

    // Most of the rows: skip straight to the full evaluation
    int ret = -1;
    if (!all && 2 * m->pending_rows.size() <= static_cast<size_t>(ng_ + 1))
      ret = eval_rows(m, need_jac);
    if (ret == -1) {
      m->n_eval_full++;
      ret = eval_full(m, need_jac);
      if (!ret) mark_rows_cached(m, need_jac, true);
    } else if (!ret) {
      mark_rows_cached(m, need_jac, false);
    }
    return ret;
  }

  int COI_CALLCONV ConoptInterface::cb_fdevalini(const double X[], const int ROWLIST[],
                                                  int MODE, int LISTSIZE, int NUMTHREAD,
                                                  int IGNERR, int* ERRCNT, int NUMVAR,
                                                  void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    const ConoptInterface& self = m->self;

    casadi_assert(NUMVAR == self.nx_, "cb_fdevalini: NUMVAR != nx_");
    casadi_assert(static_cast<size_t>(NUMVAR) == m->cached_x.size(),
                  "cb_fdevalini: NUMVAR != cached_x size");

    const bool need_jac = (MODE != 1);

    if (self.has_tape_) {
        m->cache_valid_jac.store(false, std::memory_order_relaxed);
        m->cache_valid.store(false, std::memory_order_relaxed);
        if (self.debug_) {
            casadi::uout() << "FDEvalIni (" << LISTSIZE << " rows) x:";
            for (int i = 0; i < NUMVAR; ++i)
                casadi::uout() << " " << X[i];
            casadi::uout() << "\n";
        }
        int ret = 1;
        try {
            ret = self.fdevalini_subset(m, X, ROWLIST, MODE, LISTSIZE);
        } catch (std::exception& ex) {
            casadi::uerr() << ex.what() << std::endl;
        } catch (...) {
        }
        if (ret) {
            *ERRCNT = 1;
            m->nan_encountered = true;
        } else {
            m->cache_valid_jac.store(need_jac, std::memory_order_relaxed);
            m->cache_valid.store(true, std::memory_order_release);
        }
        return 0;
    }

    // Single-row requests (typical during CONOPT's preprocessing, where rows are
    // evaluated one at a time) often repeat the same point. If so and the cache
    // already covers this MODE, skip re-evaluation. Multi-row requests always
    // re-evaluate.
    if (LISTSIZE == 1) {
        bool new_solution = !(m->stored_fun || m->stored_jac) ||
            std::memcmp(m->cached_x.data(), X, NUMVAR * sizeof(double)) != 0;
        if (new_solution) {
            m->stored_fun = false;
            m->stored_jac = false;
        } else if (m->stored_fun && (!need_jac || m->stored_jac)) {
            // FDEvalEnd cleared the validity flags; the cached data is still good.
            m->cache_valid_jac.store(m->stored_jac, std::memory_order_relaxed);
            m->cache_valid.store(true, std::memory_order_release);
            if (self.debug_)
                casadi::uout() << "FDEvalIni: point unchanged, reusing cached evaluation\n";
            return 0;
        }
    } else {
        m->stored_fun = false;
        m->stored_jac = false;
    }

    std::memcpy(m->cached_x.data(), X, NUMVAR * sizeof(double));

    if (self.debug_) {
        casadi::uout() << "FDEvalIni x:";
        for (int i = 0; i < NUMVAR; ++i)
            casadi::uout() << " " << X[i];
        casadi::uout() << "\n";
    }

    m->cache_valid_jac.store(false, std::memory_order_relaxed);
    m->cache_valid.store(false, std::memory_order_relaxed);
    try {
        int ret = self.eval_full(m, need_jac);
        if (!ret) {
            // Publish the cache only after the evaluation succeeds.
            m->cache_valid_jac.store(need_jac, std::memory_order_relaxed);
            m->cache_valid.store(true, std::memory_order_release);
            if (LISTSIZE == 1) {
                m->stored_fun = true;
                m->stored_jac = need_jac;
            }
        }
    } catch (std::exception& ex) {
        casadi::uerr() << ex.what() << std::endl;
    } catch (...) {
    }
    if (!m->cache_valid.load(std::memory_order_relaxed)) {
        *ERRCNT = 1;
        m->nan_encountered = true;
    }
    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_fd_eval(const double X[], double* G, double JAC[],
                                                int ROWNO, const int JACNUM[], int MODE,
                                                int IGNERR, int* ERRCNT, int NUMVAR,
                                                int NUMJAC, int THREAD, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    const ConoptInterface& self = m->self;

    // Acquire load: establishes happens-before with the release store in cb_fdevalini,
    // making all cache writes (including cached_jac_g and cache_valid_jac) visible here.
    if (!m->cache_valid.load(std::memory_order_acquire)) {
        *ERRCNT = 1;
        return 0;
    }

    if ((MODE == 2 || MODE == 3) &&
        !m->cache_valid_jac.load(std::memory_order_relaxed)) {
        *ERRCNT = 1;
        return 0;
    }

#ifndef NDEBUG
    casadi_assert(
      std::memcmp(X, m->cached_x.data(), NUMVAR * sizeof(double)) == 0,
      "cb_fd_eval: X does not match cached_x — CONOPT API contract violated");
#endif

    // A row that FDEvalIni did not list is evaluated here instead
    if (self.has_tape_ && ROWNO >= 0 && ROWNO <= m->ng_expanded) {
        const bool need_jac = (MODE != 1);
        int r = ROWNO == 0 ? static_cast<int>(self.ng_) : m->conopt_to_casadi[ROWNO - 1];
        if ((need_jac ? m->stamp_jac : m->stamp_fun)[r] != m->cur_stamp) {
            m->n_eval_on_demand++;
            m->pending_rows.assign(1, r);
            int ret = 1;
            try {
                ret = self.eval_rows(m, need_jac);
                if (ret == -1) {
                    m->n_eval_full++;
                    ret = self.eval_full(m, need_jac);
                    if (!ret) mark_rows_cached(m, need_jac, true);
                } else if (!ret) {
                    mark_rows_cached(m, need_jac, false);
                }
            } catch (std::exception& ex) {
                casadi::uerr() << ex.what() << std::endl;
            } catch (...) {
            }
            if (ret) {
                *ERRCNT = 1;
                return 0;
            }
        }
    }

    if (ROWNO == 0) {
        if (MODE == 1 || MODE == 3) *G = m->cached_f;
        if (MODE == 2 || MODE == 3) {
            const casadi_int* f_row = self.gradf_sp_.row();
            for (casadi_int k = 0; k < self.gradf_sp_.nnz(); ++k) {
                if (self.gradf_nlflag_[k]) {
                    JAC[f_row[k]] = m->cached_grad_f[k];
                    if (self.debug_)
                        casadi::uout() << "  df/dx[" << f_row[k] << "] = "
                                       << m->cached_grad_f[k] << "\n";
                }
            }
        }
    } else {
        casadi_assert(ROWNO >= 1 && ROWNO <= m->ng_expanded,
                      "cb_fd_eval: ROWNO out of range");
        int ci = m->conopt_to_casadi[ROWNO - 1];
        if (MODE == 1 || MODE == 3) {
            *G = m->cached_g[ci];
            if (self.debug_) {
                ConoptRowType ctype = m->conopt_type[ROWNO - 1];
                double rhs = m->conopt_rhs[ROWNO - 1];
                const char* rel = (ctype == ConoptRowType::Equality) ? "=" :
                                  (ctype == ConoptRowType::GreaterEqual) ? ">=" :
                                  (ctype == ConoptRowType::LessEqual) ? "<=" : "free";
                if (ctype == ConoptRowType::Free)
                    casadi::uout() << "  g[" << ci << "](x) = " << *G << " (free)\n";
                else
                    casadi::uout() << "  g[" << ci << "](x) = " << *G << " "
                                   << rel << " " << rhs << "\n";
            }
        }
        if (MODE == 2 || MODE == 3) {
            int base  = self.jacg_rowstart_[ci];
            int count = self.jacg_rowstart_[ci + 1] - self.jacg_rowstart_[ci];
            for (int k = 0; k < count; ++k) {
                int col = self.jacg_col_[base + k];
                double val = m->cached_jac_g[self.jacg_nzidx_[base + k]];
                JAC[col] = val;
                if (self.debug_)
                    casadi::uout() << "  dg[" << ci << "]/dx[" << col << "] = " << val << "\n";
            }
        }
    }
    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_fdevalend(int IGNERR, int* ERRCNT, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    m->cache_valid_jac.store(false, std::memory_order_relaxed);
    m->cache_valid.store(false, std::memory_order_relaxed);
    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_2dlagrstr(int HSRW[], int HSCL[], int* NODRV,
                                                  int NUMVAR, int NUMCON, int NHESS,
                                                  void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    const ConoptInterface& self = m->self;
    const casadi_int* colind = self.hesslag_sp_.colind();
    const casadi_int* row = self.hesslag_sp_.row();

    int idx = 0;
    for (int c = 0; c < NUMVAR; ++c) {
        for (casadi_int el = colind[c]; el < colind[c+1]; ++el) {
            HSRW[idx] = row[el];
            HSCL[idx] = c;
            idx++;
        }
    }
    casadi_assert(idx == NHESS,
                  "cb_2dlagrstr: idx != NHESS - Hessian nonzero count mismatch");
    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_2dlagrval(const double X[], const double U[],
                                                  const int HSRW[], const int HSCL[],
                                                  double HSVL[], int* NODRV, int NUMVAR,
                                                  int NUMCON, int NHESS, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    const ConoptInterface& self = m->self;

    // CONOPT's Lagrangian: L = SUM(r) U(r) * F(r), so d²L/dx² = SUM(r) U(r) * d²F(r)/dx².
    // CasADi computes lam_f*d²f/dx² + lam_g^T*d²g/dx², so lam_f = U[0], lam_g[ci] = U[row_ci].
    double obj_factor = U[0];

    for (int ci = 0; ci < self.ng_; ++ci) {
      int row1 = m->casadi_to_conopt_lb_row[ci];
      int row2 = m->casadi_to_conopt_ub_row[ci];
      m->hess_lam_g_[ci] = U[row1] + (row2 >= 0 ? U[row2] : 0.0);
    }

    if (self.debug_) {
        casadi::uout() << "Hessian lam_f=" << obj_factor << " lam_g:";
        for (int ci = 0; ci < self.ng_; ++ci)
            casadi::uout() << " " << m->hess_lam_g_[ci];
        casadi::uout() << "\n";
    }

    m->arg[0] = X;
    m->arg[1] = m->d_nlp.p;
    m->arg[2] = &obj_factor;
    m->arg[3] = m->hess_lam_g_.data();
    m->res[0] = HSVL;

    try {
        if (self.calc_function(m, "nlp_hess_l")) {
            *NODRV = 1;
            return 0;
        }
        if (self.debug_) {
            casadi::uout() << "Hessian values (HSVL):";
            for (int i = 0; i < NHESS; ++i)
                casadi::uout() << " " << HSVL[i];
            casadi::uout() << "\n";
        }
    } catch (std::exception& ex) {
        casadi::uerr() << "CONOPT: nlp_hess_l failed: " << ex.what() << std::endl;
        *NODRV = 1;
    } catch (...) {
        casadi::uerr() << "CONOPT: nlp_hess_l failed (unknown exception)" << std::endl;
        *NODRV = 1;
    }
    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_status(int MODSTA, int SOLSTA, int ITER,
                                               double OBJVAL, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    m->modsta = static_cast<ConoptModelStatus>(MODSTA);
    m->solsta = static_cast<ConoptSolverStatus>(SOLSTA);
    m->iter = ITER;
    m->d_nlp.objective = OBJVAL;

    const char* modsta_str;
    switch (m->modsta) {
      case ConoptModelStatus::Optimal:            modsta_str = "Optimal";                  break;
      case ConoptModelStatus::LocallyOptimal:     modsta_str = "Locally optimal";          break;
      case ConoptModelStatus::Unbounded:          modsta_str = "Unbounded";                break;
      case ConoptModelStatus::Infeasible:         modsta_str = "Infeasible";               break;
      case ConoptModelStatus::LocallyInfeasible:  modsta_str = "Locally infeasible";       break;
      case ConoptModelStatus::IntermediateInfeas: modsta_str = "Intermediate infeasible";  break;
      case ConoptModelStatus::IntermediateNonOpt: modsta_str = "Intermediate non-optimal"; break;
      case ConoptModelStatus::UnknownError:       modsta_str = "Unknown error";            break;
      case ConoptModelStatus::ErrorNoSolution:    modsta_str = "Error: no solution";       break;
      default:                                    modsta_str = "Unknown model status";     break;
    }

    const char* solsta_str;
    switch (m->solsta) {
      case ConoptSolverStatus::NormalCompletion:
        solsta_str = "Normal completion";                   break;
      case ConoptSolverStatus::IterationLimit:
        solsta_str = "Iteration limit";                     break;
      case ConoptSolverStatus::TimeLimit:
        solsta_str = "Time limit";                          break;
      case ConoptSolverStatus::TerminatedBySolver:
        solsta_str = "Terminated by solver";                break;
      case ConoptSolverStatus::EvalErrorLimit:
        solsta_str = "Evaluation error limit";              break;
      case ConoptSolverStatus::UserInterrupt:
        solsta_str = "User interrupt";                      break;
      case ConoptSolverStatus::SetupFailure:
        solsta_str = "Setup failure";                       break;
      case ConoptSolverStatus::MajorSolverError:
        solsta_str = "Major solver error";                  break;
      case ConoptSolverStatus::MajorSolverErrorFeas:
        solsta_str = "Major solver error (feasible point)"; break;
      case ConoptSolverStatus::SystemError:
        solsta_str = "System error";                        break;
      case ConoptSolverStatus::QuickModeTermination:
        solsta_str = "Quick Mode termination";              break;
      default:
        solsta_str = "Unknown solver status";               break;
    }

    m->return_status = std::string(modsta_str) + " / " + std::string(solsta_str);
    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_solution(const double XVAL[], const double XMAR[],
                                                 const int XBAS[], const int XSTA[],
                                                 const double YVAL[], const double YMAR[],
                                                 const int YBAS[], const int YSTA[],
                                                 int NUMVAR, int NUMCON, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);

    casadi_copy(XVAL, NUMVAR, m->d_nlp.z);

    if (m->self.debug_) {
        casadi::uout() << "Solution x:";
        for (int i = 0; i < NUMVAR; ++i)
            casadi::uout() << " " << XVAL[i];
        casadi::uout() << "\n";
    }

    // Use the first row of each constraint and restore constants absorbed into the RHS.
    for (casadi_int ci = 0; ci < m->self.ng_; ++ci)
      m->d_nlp.z[NUMVAR + ci] = YVAL[m->casadi_to_conopt_lb_row[ci]] + m->row_const_[ci];

    // Variable marginals: CONOPT shadow prices = -CasADi lam_x
    for (int i = 0; i < NUMVAR; ++i)
      m->d_nlp.lam[i] = -XMAR[i];

    // Constraint marginals: for range constraints sum both rows' shadow prices
    // (only the active bound has a non-zero YMAR; summing is always safe)
    for (casadi_int ci = 0; ci < m->self.ng_; ++ci) {
      int row1 = m->casadi_to_conopt_lb_row[ci];
      int row2 = m->casadi_to_conopt_ub_row[ci];
      double ymar = YMAR[row1] + (row2 >= 0 ? YMAR[row2] : 0.0);
      m->d_nlp.lam[NUMVAR + ci] = -ymar;
    }

    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_message(int SMSG, int DMSG, int NMSG, char* MSGV[],
                                                void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);

    int message_length = SMSG;
    if (m->self.debug_) message_length = std::max(message_length, std::max(DMSG, NMSG));

    for (int i = 0; i < message_length; ++i) {
        if (MSGV[i] != nullptr) {
            casadi::uout() << MSGV[i] << std::endl;
        }
    }
    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_errmsg(int ROWNO, int COLNO, int POSNO,
                                               const char* MSG, void* USRMEM) {
    if (MSG == nullptr) return 0;

    std::string prefix = "CONOPT Error: ";
    if (COLNO == -1 && ROWNO >= 0) {
        prefix += "Row " + std::to_string(ROWNO) + " - ";
    } else if (ROWNO == -1 && COLNO >= 0) {
        prefix += "Column " + std::to_string(COLNO) + " - ";
    } else if (ROWNO >= 0 && COLNO >= 0) {
        if (POSNO >= 0) {
            prefix += "Jacobian Pos " + std::to_string(POSNO) + " (Row " +
                      std::to_string(ROWNO) + ", Col " + std::to_string(COLNO) + ") - ";
        } else if (POSNO == -1) {
            prefix += "Pair (Row " + std::to_string(ROWNO) + ", Col " +
                      std::to_string(COLNO) + ") - ";
        }
    }

    casadi::uerr() << prefix << MSG << std::endl;
    return 0;
  }
}  // namespace casadi
