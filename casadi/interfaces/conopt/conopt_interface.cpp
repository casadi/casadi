#include "conopt_interface.hpp"
#include "casadi/core/casadi_misc.hpp"
#include "casadi/core/casadi_interrupt.hpp"
#include <cassert>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <limits>

namespace casadi {

  extern "C" int CASADI_NLPSOL_CONOPT_EXPORT casadi_register_nlpsol_conopt(Nlpsol::Plugin* plugin) {
    plugin->creator = ConoptInterface::creator;
    plugin->name = "conopt";
    plugin->doc = "CONOPT Interface";
    plugin->version = CASADI_VERSION;
    plugin->options = &ConoptInterface::options_;
    plugin->deserialize = &ConoptInterface::deserialize;
    return 0;
  }

  extern "C" void CASADI_NLPSOL_CONOPT_EXPORT casadi_load_nlpsol_conopt() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_conopt);
  }

  const Options ConoptInterface::options_ = {{&Nlpsol::options_}, {
      {"exact_hessian", {OT_BOOL, "Provide exact Hessian to CONOPT"}},
      {"warm_start", {OT_BOOL, "Warm-start CONOPT using multipliers from a prior solve to infer basis status (IniStat=2)"}},
      {"conopt", {OT_DICT, "Options to be passed to CONOPT"}},
      {"optfile", {OT_STRING, "Path to a CONOPT option file (for string-valued CR-cells such as Algorithm)"}},
      {"debug", {OT_BOOL, "Print debug output: constraint values at each FDEval, solution vector, and option echo"}}
  }};

  ConoptInterface::ConoptInterface(const std::string& name, const Function& nlp)
      : Nlpsol(name, nlp) {}

  ConoptInterface::~ConoptInterface() { clear_mem(); }

  void ConoptInterface::init(const Dict& opts) {
    Nlpsol::init(opts);

    // Extract native options
    warm_start_ = false;
    debug_ = false;
    for (auto&& op : opts) {
      if (op.first == "conopt") opts_ = op.second;
      else if (op.first == "optfile") optfile_ = op.second.to_string();
      else if (op.first == "warm_start") warm_start_ = op.second.to_bool();
      else if (op.first == "debug") debug_ = op.second.to_bool();
    }

    create_function("nlp_f", {"x", "p"}, {"f"});
    create_function("nlp_g", {"x", "p"}, {"g"});
    Function gradf_fcn = create_function("nlp_grad_f", {"x", "p"}, {"f", "grad:f:x"});
    gradf_sp_ = gradf_fcn.sparsity_out(1);

    Function jacg_fcn = create_function("nlp_jac_g", {"x", "p"}, {"g", "jac:g:x"});
    jacg_sp_ = jacg_fcn.sparsity_out(1);

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

    const casadi_int* g_colind = jacg_sp_.colind();
    const casadi_int* g_row = jacg_sp_.row();

    // Build row-indexed CSR structure over nonlinear entries for fast Jacobian scatter in cb_fd_eval
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
        int r = (int)g_row[el];
        int pos = jacg_rowstart_[r] + fill_pos[r]++;
        jacg_nzidx_[pos] = (int)el;
        jacg_col_[pos]   = c;
      }
    }
  }

  // --- Serialization & Deserialization --- //
  ConoptInterface::ConoptInterface(DeserializingStream& s) : Nlpsol(s) {
    s.version("ConoptInterface", 1);
    s.unpack("ConoptInterface::exact_hessian", exact_hessian_);
    s.unpack("ConoptInterface::opts", opts_);
    s.unpack("ConoptInterface::gradf_sp", gradf_sp_);
    s.unpack("ConoptInterface::jacg_sp", jacg_sp_);
    s.unpack("ConoptInterface::hesslag_sp", hesslag_sp_);
    s.unpack("ConoptInterface::optfile", optfile_);
    s.unpack("ConoptInterface::warm_start", warm_start_);
    s.unpack("ConoptInterface::debug", debug_);

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
        int r = (int)g_row[el];
        int pos = jacg_rowstart_[r] + fill_pos[r]++;
        jacg_nzidx_[pos] = (int)el;
        jacg_col_[pos]   = c;
      }
  }

  void ConoptInterface::serialize_body(SerializingStream &s) const {
    Nlpsol::serialize_body(s);
    s.version("ConoptInterface", 1);
    s.pack("ConoptInterface::exact_hessian", exact_hessian_);
    s.pack("ConoptInterface::opts", opts_);
    s.pack("ConoptInterface::gradf_sp", gradf_sp_);
    s.pack("ConoptInterface::jacg_sp", jacg_sp_);
    s.pack("ConoptInterface::hesslag_sp", hesslag_sp_);
    s.pack("ConoptInterface::optfile", optfile_);
    s.pack("ConoptInterface::warm_start", warm_start_);
    s.pack("ConoptInterface::debug", debug_);
    // Derived arrays (gradf_col_flag_, CSR) are rebuilt on deserialization
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
    m->row_nnz.assign(ng_, 0);
    // Most problems have few or no range constraints, so reserve a fraction of
    // the worst case (all rows range); solve() grows these
    // vectors on demand as rows are actually expanded.
    casadi_int initial_row_reserve = std::max<casadi_int>(1, ng_ / 4);
    m->conopt_to_casadi.reserve(initial_row_reserve);
    m->conopt_type.reserve(initial_row_reserve);
    m->conopt_rhs.reserve(initial_row_reserve);
    if (has_linear_jac_)
      m->const_jac_vals.resize(jacg_sp_.nnz(), 0.0);
    if (has_linear_gradf_)
      m->gradf_const_vals.resize(gradf_sp_.nnz(), 0.0);

    if (COI_Create(&m->cntvect) != 0 || m->cntvect == nullptr) {
      casadi::uerr() << "CONOPT: COI_Create failed" << std::endl;
      return 1;
    }

    if (warm_start_) COIDEF_IniStat(m->cntvect, 2);

    COIDEF_NumVar(m->cntvect, nx_);
    // NumCon, NumNz, NumNlNz are set in solve() because range-constraint expansion
    // can change them between calls.

    COIDEF_ObjCon(m->cntvect, 0);
    COIDEF_OptDir(m->cntvect, -1);

    // Handle Options
    m->custom_options.clear();
    for (auto&& op : opts_) {
        // Explictly catch C API options defined in conopt.h
        if (op.first == "itlim") COIDEF_ItLim(m->cntvect, op.second.to_int());
        else if (op.first == "errlim") COIDEF_ErrLim(m->cntvect, op.second.to_int());
        else if (op.first == "reslim") COIDEF_ResLim(m->cntvect, op.second.to_double());
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
        COIDEF_2DLagrSize(m->cntvect, &ConoptInterface::cb_2dlagrsize);
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

  void ConoptInterface::set_work(void* mem, const double**& arg, double**& res, casadi_int*& iw, double*& w) const {
    Nlpsol::set_work(mem, arg, res, iw, w);
  }

  int ConoptInterface::solve(void* mem) const {
    auto m = static_cast<ConoptMemory*>(mem);
    m->cache_valid     = false;
    m->cache_valid_jac = false;
    m->cached_f        = 0.0;
    m->nan_encountered = false;

    // Build the per-solve constraint expansion (splits range constraints into two rows)
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

    // conopt_to_casadi/conopt_type/conopt_rhs always grow in lockstep, so a single
    // capacity check (on conopt_to_casadi) is enough to decide whether to grow all
    // three. Growth is by 0.25 of the CasADi rows not yet processed, rather than
    // jumping straight to the worst-case (all-range) size.
    auto ensure_row_capacity = [m](casadi_int remaining_rows) {
      if (m->conopt_to_casadi.size() == m->conopt_to_casadi.capacity()) {
        casadi_int new_cap = m->conopt_to_casadi.capacity() +
                              std::max<casadi_int>(std::min<casadi_int>(remaining_rows, 10), remaining_rows / 4);
        m->conopt_to_casadi.reserve(new_cap);
        m->conopt_type.reserve(new_cap);
        m->conopt_rhs.reserve(new_cap);
      }
    };

    for (casadi_int i = 0; i < ng_; ++i) {
      double lbg = m->d_nlp.lbz[nx_ + i];
      double ubg = m->d_nlp.ubz[nx_ + i];
      bool is_range = !std::isinf(lbg) && !std::isinf(ubg) && lbg != ubg;

      // CONOPT row 0 is reserved for the objective, so constraint rows start at 1
      // (arrays are still plain 0-based C arrays; only the row *numbering* is offset).
      m->casadi_to_conopt_lb_row[i] = (int)(ng_expanded + 1);
      ensure_row_capacity(ng_ - i);
      m->conopt_to_casadi.push_back((int)i);
      if (lbg == ubg) {
        m->conopt_type.push_back(0);  m->conopt_rhs.push_back(lbg);
      } else if (!std::isinf(lbg) && std::isinf(ubg)) {
        m->conopt_type.push_back(1);  m->conopt_rhs.push_back(lbg);
      } else if (std::isinf(lbg) && !std::isinf(ubg)) {
        m->conopt_type.push_back(2);  m->conopt_rhs.push_back(ubg);
      } else if (std::isinf(lbg) && std::isinf(ubg)) {
        m->conopt_type.push_back(3);  m->conopt_rhs.push_back(0.0);
      } else {
        m->conopt_type.push_back(1);  m->conopt_rhs.push_back(lbg);  // range: >= row
      }
      ng_expanded++;

      if (is_range) {
        m->casadi_to_conopt_ub_row[i] = (int)(ng_expanded + 1);
        ensure_row_capacity(ng_ - i);
        m->conopt_to_casadi.push_back((int)i);
        m->conopt_type.push_back(2);  m->conopt_rhs.push_back(ubg);  // <= row
        ng_expanded++;
        numnz += m->row_nnz[i];
      }
    }
    casadi_assert(ng_expanded <= std::numeric_limits<int>::max(), "ng_expanded overflows int");
    m->ng_expanded = (int)ng_expanded;

    // Evaluate Jacobian at the initial point to obtain values for constant entries
    if (has_linear_jac_) {
      m->arg[0] = m->d_nlp.z;
      m->arg[1] = m->d_nlp.p;
      m->res[0] = m->cached_g.data();
      m->res[1] = m->const_jac_vals.data();
      try {
        calc_function(m, "nlp_jac_g");
      } catch (std::exception& ex) {
        casadi::uerr() << "CONOPT: initial evaluation failed: " << ex.what() << std::endl;
        return 1;
      } catch (...) {
        casadi::uerr() << "CONOPT: initial evaluation failed (unknown exception)" << std::endl;
        return 1;
      }
    }

    // Adjust conopt_rhs for fully linear rows to absorb constant terms.
    // CONOPT's pre-triangular preprocessor uses only VALUE[] + RHS and never calls
    // FDEval during that phase, so any constant in a fully linear constraint
    // (e.g. x[9] - 3*x[6] + 133 = 0 → RHS must be -133, not 0) must be moved
    // to the RHS here. For mixed (linear+nonlinear) rows, CONOPT calls FDEval and
    // gets the full function value including the constant, so no adjustment needed.
    if (has_linear_jac_) {
      const casadi_int* g_colind_c = jacg_sp_.colind();
      const casadi_int* g_row_c    = jacg_sp_.row();

      // Accumulate the linear part of G at x0 per row: sum_j a_j * x0_j
      std::vector<double> linear_at_x0(ng_, 0.0);
      for (int c = 0; c < nx_; ++c) {
        for (casadi_int el = g_colind_c[c]; el < g_colind_c[c + 1]; ++el) {
          if (jacg_nlflag_[el] == 0)
            linear_at_x0[g_row_c[el]] += m->const_jac_vals[el] * m->d_nlp.z[c];
        }
      }

      for (int ci = 0; ci < ng_; ++ci) {
        // Only adjust fully linear rows (no nonlinear Jacobian entries)
        if (jacg_rowstart_[ci + 1] != jacg_rowstart_[ci]) continue;
        double constant = m->cached_g[ci] - linear_at_x0[ci];
        if (std::abs(constant) < 1e-14) continue;
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
        calc_function(m, "nlp_grad_f");
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
          int ci = (int)g_row_s[el];
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
    m->numnz_expanded = (int)numnz;

    COIDEF_NumCon(m->cntvect, (int)(ng_expanded + 1));
    COIDEF_NumNz(m->cntvect, (int)numnz);
    COIDEF_NumNlNz(m->cntvect, (int)num_nl_nz);

    int ret = COI_Solve(m->cntvect);

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
    return stats;
  }

  // --- Dynamic Option Callback --- //
  int COI_CALLCONV ConoptInterface::cb_option(int NCALL, double* RVAL, int* IVAL, int* LVAL, char* NAME, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);

    if (NCALL >= (int)m->custom_options.size()) {
        NAME[0] = '\0';
        return 0;
    }

    auto& opt = m->custom_options[NCALL];
    std::strcpy(NAME, opt.first.c_str());

    if (opt.second.is_double()) {
        *RVAL = opt.second.to_double();
        if (m->self.debug_) casadi::uout() << "CONOPT option: " << opt.first << " = " << *RVAL << std::endl;
    } else if (opt.second.is_int()) {
        *IVAL = opt.second.to_int();
        if (m->self.debug_) casadi::uout() << "CONOPT option: " << opt.first << " = " << *IVAL << std::endl;
    } else if (opt.second.is_bool()) {
        *LVAL = opt.second.to_bool() ? 1 : 0;
        if (m->self.debug_) casadi::uout() << "CONOPT option: " << opt.first << " = " << (opt.second.to_bool() ? "true" : "false") << std::endl;
    }
    else if (opt.second.is_string()) {
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
  int COI_CALLCONV ConoptInterface::cb_progress(int LEN_INT, const int INTX[], int LEN_RL, const double RL[], const double X[], void* USRMEM) {
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

  int COI_CALLCONV ConoptInterface::cb_read_matrix(double LOWER[], double CURR[], double UPPER[], int VSTA[], int TYPEX[], double RHS[], int ESTA[], int COLSTA[], int ROWNO[], double VALUE[], int NLFLAG[], int NUMVAR, int NUMCON, int NUMNZ, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    const ConoptInterface& self = m->self;

    // Variable bounds and initial point (clamped to [lb, ub])
    for (int i = 0; i < NUMVAR; ++i) {
      double lb = m->d_nlp.lbz[i];
      double ub = m->d_nlp.ubz[i];
      if (!std::isinf(lb)) LOWER[i] = lb;
      if (!std::isinf(ub)) UPPER[i] = ub;
      double x0 = m->d_nlp.z[i];
      if (!std::isinf(lb)) x0 = std::max(x0, lb);
      if (!std::isinf(ub)) x0 = std::min(x0, ub);
      CURR[i] = x0;
    }

    // Constraint types and RHS (row 0 = objective, rows 1..ng_expanded = constraints)
    TYPEX[0] = 3;
    RHS[0]   = 0.0;
    for (int r = 0; r < m->ng_expanded; ++r) {
      TYPEX[r + 1] = m->conopt_type[r];
      RHS[r + 1]   = m->conopt_rhs[r];
    }

    if (self.warm_start_) {
      // Always zero-initialise so that skipping the fill below is safe.
      std::fill(VSTA, VSTA + NUMVAR, 0);
      std::fill(ESTA, ESTA + NUMCON, 0);
      ESTA[0] = 3;  // objective row always superbasic

      const double* lam = m->d_nlp.lam;
      bool all_zero = std::all_of(lam, lam + self.nx_ + self.ng_,
                                  [](double v) { return v == 0.0; });
      if (!all_zero) {
        for (int i = 0; i < NUMVAR; ++i) {
          double lbi = m->d_nlp.lbz[i];
          double ubi = m->d_nlp.ubz[i];
          double xi  = m->d_nlp.z[i];
          double li  = lam[i];
          if (!std::isinf(lbi) && std::fabs(xi - lbi) < 1e-8 && li <= 0.0)
            VSTA[i] = 0;
          else if (!std::isinf(ubi) && std::fabs(xi - ubi) < 1e-8 && li >= 0.0)
            VSTA[i] = 1;
          else
            VSTA[i] = 2;
        }

        for (int r = 0; r < m->ng_expanded; ++r) {
          int ci        = m->conopt_to_casadi[r];
          double lam_ci = lam[NUMVAR + ci];
          int row1      = m->casadi_to_conopt_lb_row[ci];
          int row2      = m->casadi_to_conopt_ub_row[ci];
          if (r + 1 == row1) {
            int ctype = m->conopt_type[r];  // type of this CONOPT expanded row
            if (ctype == 0) {               // equality: both sides, just mark basic
              ESTA[r + 1] = 2;
            } else if (ctype == 1) {        // >= row: active when lam_ci < 0
              ESTA[r + 1] = (lam_ci < 0.0) ? 0 : 2;
            } else if (ctype == 2) {        // <= row (pure <= stored as lb_row): active when lam_ci > 0
              ESTA[r + 1] = (lam_ci > 0.0) ? 1 : 2;
            } else {                        // free row (type=3): superbasic
              ESTA[r + 1] = 3;
            }
          } else if (r + 1 == row2) {
            // range ub side is active when lam_ci > 0
            ESTA[r + 1] = (lam_ci > 0.0) ? 1 : 2;
          } else {
            ESTA[r + 1] = 2;
          }
        }
      }
    }

    // Jacobian structure — built live from jacg_sp_ with range-row duplication
    const casadi_int* g_colind = self.jacg_sp_.colind();
    const casadi_int* g_row    = self.jacg_sp_.row();
    int nz = 0;
    for (int c = 0; c < NUMVAR; ++c) {
      COLSTA[c] = nz;
      if (self.gradf_col_flag_[c]) {
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
      for (casadi_int el = g_colind[c]; el < g_colind[c+1]; ++el) {
        int ci = (int)g_row[el];
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

    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_fdevalini(const double X[], const int ROWLIST[], int MODE, int LISTSIZE, int NUMTHREAD, int IGNERR, int* ERRCNT, int NUMVAR, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    const ConoptInterface& self = m->self;

    std::memcpy(m->cached_x.data(), X, NUMVAR * sizeof(double));

    if (self.debug_) {
        casadi::uout() << "FDEvalIni x:";
        for (int i = 0; i < NUMVAR; ++i)
            casadi::uout() << " " << X[i];
        casadi::uout() << "\n";
    }

    const bool need_jac = (MODE != 1);

    try {
        m->arg[0] = m->cached_x.data();
        m->arg[1] = m->d_nlp.p;

        m->res[0] = &m->cached_f;
        m->res[1] = m->cached_grad_f.data();
        self.calc_function(m, "nlp_grad_f");

        // When MODE==1, only function values (G) are needed — skip Jacobian computation.
        m->res[0] = m->cached_g.data();
        m->res[1] = need_jac ? m->cached_jac_g.data() : nullptr;
        self.calc_function(m, "nlp_jac_g");

        // cache_valid_jac is stored before cache_valid (release).  The release on
        // cache_valid orders both stores, so readers need only an acquire on cache_valid.
        m->cache_valid_jac.store(need_jac, std::memory_order_relaxed);
        m->cache_valid.store(true, std::memory_order_release);
    } catch (std::exception& ex) {
        casadi::uerr() << ex.what() << std::endl;
        *ERRCNT = 1;
        m->nan_encountered = true;
        m->cache_valid_jac.store(false, std::memory_order_relaxed);
        m->cache_valid.store(false, std::memory_order_relaxed);
    } catch (...) {
        *ERRCNT = 1;
        m->nan_encountered = true;
        m->cache_valid_jac.store(false, std::memory_order_relaxed);
        m->cache_valid.store(false, std::memory_order_relaxed);
    }

    // Detect silent NaN from CasADi (e.g. sqrt of negative number)
    if (m->cache_valid.load(std::memory_order_relaxed)) {
        bool has_nan = std::isnan(m->cached_f);
        if (!has_nan) {
            for (double v : m->cached_g) if (std::isnan(v)) { has_nan = true; break; }
        }
        if (!has_nan) {
            for (double v : m->cached_grad_f) if (std::isnan(v)) { has_nan = true; break; }
        }
        if (!has_nan && need_jac) {
            for (double v : m->cached_jac_g) if (std::isnan(v)) { has_nan = true; break; }
        }
        if (has_nan) {
            *ERRCNT = 1;
            m->nan_encountered = true;
            m->cache_valid_jac.store(false, std::memory_order_relaxed);
            m->cache_valid.store(false, std::memory_order_relaxed);
        }
    }
    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_fd_eval(const double X[], double* G, double JAC[], int ROWNO, const int JACNUM[], int MODE, int IGNERR, int* ERRCNT, int NUMVAR, int NUMJAC, int THREAD, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    const ConoptInterface& self = m->self;

    // Acquire load: establishes happens-before with the release store in cb_fdevalini,
    // making all cache writes (including cached_jac_g and cache_valid_jac) visible here.
    if (!m->cache_valid.load(std::memory_order_acquire)) {
        *ERRCNT = 1;
        return 0;
    }

#ifndef NDEBUG
    casadi_assert(
      std::memcmp(X, m->cached_x.data(), NUMVAR * sizeof(double)) == 0,
      "cb_fd_eval: X does not match cached_x — CONOPT API contract violated");
#endif

    if (ROWNO == 0) {
        if (MODE == 1 || MODE == 3) *G = m->cached_f;
        if (MODE == 2 || MODE == 3) {
            const casadi_int* f_row = self.gradf_sp_.row();
            for (casadi_int k = 0; k < self.gradf_sp_.nnz(); ++k) {
                if (self.gradf_nlflag_[k]) {
                    JAC[f_row[k]] = m->cached_grad_f[k];
                    if (self.debug_)
                        casadi::uout() << "  df/dx[" << f_row[k] << "] = " << m->cached_grad_f[k] << "\n";
                }
            }
        }
    } else {
        int ci = m->conopt_to_casadi[ROWNO - 1];
        if (MODE == 1 || MODE == 3) {
            *G = m->cached_g[ci];
            if (self.debug_) {
                int ctype = m->conopt_type[ROWNO - 1];
                double rhs = m->conopt_rhs[ROWNO - 1];
                const char* rel = (ctype == 0) ? "=" : (ctype == 1) ? ">=" : (ctype == 2) ? "<=" : "free";
                if (ctype == 3)
                    casadi::uout() << "  g[" << ci << "](x) = " << *G << " (free)\n";
                else
                    casadi::uout() << "  g[" << ci << "](x) = " << *G << " " << rel << " " << rhs << "\n";
            }
        }
        if (MODE == 2 || MODE == 3) {
            // Guard against stale Jacobian: cache_valid_jac is false when cb_fdevalini
            // was called with MODE==1 and skipped Jacobian computation.
            if (!m->cache_valid_jac.load(std::memory_order_relaxed)) {
                *ERRCNT = 1;
                return 0;
            }
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

  int COI_CALLCONV ConoptInterface::cb_2dlagrsize(int* NODRV, int NUMVAR, int NUMCON, int* NHESS, int MAXHESS, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    const ConoptInterface& self = m->self;
    *NHESS = self.hesslag_sp_.nnz();
    if (*NHESS > MAXHESS) {
        *NODRV = 1;
    }
    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_2dlagrstr(int HSRW[], int HSCL[], int* NODRV, int NUMVAR, int NUMCON, int NHESS, void* USRMEM) {
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
    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_2dlagrval(const double X[], const double U[], const int HSRW[], const int HSCL[], double HSVL[], int* NODRV, int NUMVAR, int NUMCON, int NHESS, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    const ConoptInterface& self = m->self;

    // CONOPT's Lagrangian: L = SUM(r) U(r) * F(r), so d²L/dx² = SUM(r) U(r) * d²F(r)/dx².
    // CasADi computes lam_f*d²f/dx² + lam_g^T*d²g/dx², so lam_f = U[0], lam_g[ci] = U[row_ci].
    // No sign flip: CONOPT's U is the Lagrangian weight directly, not the shadow price (YMAR).
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
        self.calc_function(m, "nlp_hess_l");
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

  int COI_CALLCONV ConoptInterface::cb_status(int MODSTA, int SOLSTA, int ITER, double OBJVAL, void* USRMEM) {
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
      case ConoptSolverStatus::NormalCompletion:     solsta_str = "Normal completion";                   break;
      case ConoptSolverStatus::IterationLimit:       solsta_str = "Iteration limit";                     break;
      case ConoptSolverStatus::TimeLimit:            solsta_str = "Time limit";                          break;
      case ConoptSolverStatus::TerminatedBySolver:   solsta_str = "Terminated by solver";                break;
      case ConoptSolverStatus::EvalErrorLimit:       solsta_str = "Evaluation error limit";              break;
      case ConoptSolverStatus::UserInterrupt:        solsta_str = "User interrupt";                      break;
      case ConoptSolverStatus::SetupFailure:         solsta_str = "Setup failure";                       break;
      case ConoptSolverStatus::MajorSolverError:     solsta_str = "Major solver error";                  break;
      case ConoptSolverStatus::MajorSolverErrorFeas: solsta_str = "Major solver error (feasible point)"; break;
      case ConoptSolverStatus::SystemError:          solsta_str = "System error";                        break;
      case ConoptSolverStatus::QuickModeTermination: solsta_str = "Quick Mode termination";              break;
      default:                                       solsta_str = "Unknown solver status";               break;
    }

    m->return_status = std::string(modsta_str) + " / " + std::string(solsta_str);
    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_solution(const double XVAL[], const double XMAR[], const int XBAS[], const int XSTA[], const double YVAL[], const double YMAR[], const int YBAS[], const int YSTA[], int NUMVAR, int NUMCON, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);

    casadi_copy(XVAL, NUMVAR, m->d_nlp.z);

    if (m->self.debug_) {
        casadi::uout() << "Solution x:";
        for (int i = 0; i < NUMVAR; ++i)
            casadi::uout() << " " << XVAL[i];
        casadi::uout() << "\n";
    }

    // Constraint values: use the first CONOPT row for each CasADi constraint
    // (both rows carry the same function value for range constraints)
    for (casadi_int ci = 0; ci < m->self.ng_; ++ci)
      m->d_nlp.z[NUMVAR + ci] = YVAL[m->casadi_to_conopt_lb_row[ci]];

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

  int COI_CALLCONV ConoptInterface::cb_message(int SMSG, int DMSG, int NMSG, char* MSGV[], void* USRMEM) {
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

  int COI_CALLCONV ConoptInterface::cb_errmsg(int ROWNO, int COLNO, int POSNO, const char* MSG, void* USRMEM) {
    if (MSG == nullptr) return 0;

    std::string prefix = "CONOPT Error: ";
    if (COLNO == -1 && ROWNO >= 0) {
        prefix += "Row " + std::to_string(ROWNO) + " - ";
    } else if (ROWNO == -1 && COLNO >= 0) {
        prefix += "Column " + std::to_string(COLNO) + " - ";
    } else if (ROWNO >= 0 && COLNO >= 0) {
        if (POSNO >= 0) {
            prefix += "Jacobian Pos " + std::to_string(POSNO) + " (Row " + std::to_string(ROWNO) + ", Col " + std::to_string(COLNO) + ") - ";
        } else if (POSNO == -1) {
            prefix += "Pair (Row " + std::to_string(ROWNO) + ", Col " + std::to_string(COLNO) + ") - ";
        }
    }

    casadi::uerr() << prefix << MSG << std::endl;
    return 0;
  }
}
