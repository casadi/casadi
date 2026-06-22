#include "conopt_interface.hpp"
#include "casadi/core/casadi_misc.hpp"
#include "casadi/core/casadi_interrupt.hpp"
#include <cmath>
#include <cstring>
#include <algorithm>

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
      {"conopt", {OT_DICT, "Options to be passed to CONOPT"}}
  }};

  ConoptInterface::ConoptInterface(const std::string& name, const Function& nlp)
      : Nlpsol(name, nlp) {}

  ConoptInterface::~ConoptInterface() { clear_mem(); }

  void ConoptInterface::init(const Dict& opts) {
    Nlpsol::init(opts);

    // Extract native options
    for (auto&& op : opts) {
      if (op.first == "conopt") opts_ = op.second;
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
      alloc_w(hesslag_sp_.nnz(), true);
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
      for (casadi_int el = 0; el < dgradf_dx.nnz(); ++el)
        gradf_nlflag_[dg_row[el]] = 1;
      has_linear_gradf_ = std::any_of(gradf_nlflag_.begin(), gradf_nlflag_.end(),
                                       [](int f) { return f == 0; });
    }

    const casadi_int* g_colind = jacg_sp_.colind();
    const casadi_int* g_row = jacg_sp_.row();

    // Build row-indexed structure over nonlinear entries only (aligns with CONOPT's JACNUM)
    casadi_int nnz_g = jacg_sp_.nnz();
    jacg_rowstart_.assign(ng_ + 1, 0);
    for (casadi_int el = 0; el < nnz_g; ++el)
      if (jacg_nlflag_[el]) jacg_rowstart_[g_row[el] + 1]++;
    for (int r = 0; r < ng_; ++r)
      jacg_rowstart_[r + 1] += jacg_rowstart_[r];
    jacg_nzidx_.resize(jacg_rowstart_[ng_]);
    std::vector<int> fill_pos(ng_, 0);
    for (int c = 0; c < nx_; ++c) {
      for (casadi_int el = g_colind[c]; el < g_colind[c+1]; ++el) {
        if (!jacg_nlflag_[el]) continue;
        int r = (int)g_row[el];
        jacg_nzidx_[jacg_rowstart_[r] + fill_pos[r]++] = (int)el;
      }
    }

    alloc_w(1, true);
    alloc_w(gradf_sp_.nnz(), true);
    alloc_w(ng_, true);
    alloc_w(jacg_sp_.nnz(), true);
  }

  // --- Serialization & Deserialization --- //
  ConoptInterface::ConoptInterface(DeserializingStream& s) : Nlpsol(s) {
    s.version("ConoptInterface", 1);
    s.unpack("ConoptInterface::exact_hessian", exact_hessian_);
    s.unpack("ConoptInterface::opts", opts_);
    s.unpack("ConoptInterface::gradf_sp", gradf_sp_);
    s.unpack("ConoptInterface::jacg_sp", jacg_sp_);
    s.unpack("ConoptInterface::hesslag_sp", hesslag_sp_);

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
    std::vector<int> fill_pos(ng_, 0);
    for (int c = 0; c < nx_; ++c)
      for (casadi_int el = g_colind[c]; el < g_colind[c+1]; ++el) {
        if (!jacg_nlflag_[el]) continue;
        int r = (int)g_row[el];
        jacg_nzidx_[jacg_rowstart_[r] + fill_pos[r]++] = (int)el;
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
    // Derived arrays (gradf_col_flag_, CSR) are rebuilt on deserialization
  }

  ConoptMemory::ConoptMemory(const ConoptInterface& interface)
      : self(interface), NlpsolMemory(), cntvect(nullptr), modsta(0), solsta(0),
        iter(0), return_status("Unset"),
        cache_valid(false), nan_encountered(false),
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
    if (has_linear_jac_)
      m->const_jac_vals.resize(jacg_sp_.nnz(), 0.0);
    if (has_linear_gradf_)
      m->gradf_const_vals.resize(gradf_sp_.nnz(), 0.0);

    COI_Create(&m->cntvect);

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
        else {
            // Push to custom options for the cb_option callback to feed dynamically
            m->custom_options.push_back(op);
        }
    }
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
    m->cache_valid = false;
    m->nan_encountered = false;

    // Build the per-solve constraint expansion (splits range constraints into two rows)
    m->conopt_to_casadi.clear();
    m->casadi_to_conopt_ub_row.assign(ng_, -1);
    m->conopt_type.clear();
    m->conopt_rhs.clear();

    // Compute total nnz per CasADi row (needed for range-constraint NZ duplication)
    std::vector<int> row_nnz(ng_, 0);
    {
      const casadi_int* g_row_s = jacg_sp_.row();
      for (casadi_int el = 0; el < (casadi_int)jacg_sp_.nnz(); ++el)
        row_nnz[g_row_s[el]]++;
    }

    int ng_expanded = 0;
    // Base NZ: objective gradient columns + all constraint Jacobian entries
    int numnz = (int)gradf_sp_.nnz() + (int)jacg_sp_.nnz();

    for (casadi_int i = 0; i < ng_; ++i) {
      double lbg = m->d_nlp.lbz[nx_ + i];
      double ubg = m->d_nlp.ubz[nx_ + i];
      bool is_range = !std::isinf(lbg) && !std::isinf(ubg) && lbg != ubg;

      m->casadi_to_conopt_lb_row[i] = ng_expanded + 1;  // 1-indexed CONOPT row
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
        m->casadi_to_conopt_ub_row[i] = ng_expanded + 1;
        m->conopt_to_casadi.push_back((int)i);
        m->conopt_type.push_back(2);  m->conopt_rhs.push_back(ubg);  // <= row
        ng_expanded++;
        numnz += row_nnz[i];
      }
    }
    m->ng_expanded    = ng_expanded;
    m->numnz_expanded = numnz;

    // Evaluate Jacobian at the initial point to obtain values for constant entries
    if (has_linear_jac_) {
      m->arg[0] = m->d_nlp.z;
      m->arg[1] = m->d_nlp.p;
      m->res[0] = m->cached_g.data();
      m->res[1] = m->const_jac_vals.data();
      try {
        calc_function(m, "nlp_jac_g");
      } catch (...) {
        std::fill(m->const_jac_vals.begin(), m->const_jac_vals.end(), 0.0);
      }
    }

    // Evaluate objective gradient at initial point for constant entries
    if (has_linear_gradf_) {
      m->arg[0] = m->d_nlp.z;
      m->arg[1] = m->d_nlp.p;
      m->res[0] = &m->cached_f;
      m->res[1] = m->gradf_const_vals.data();
      try {
        calc_function(m, "nlp_grad_f");
      } catch (...) {
        std::fill(m->gradf_const_vals.begin(), m->gradf_const_vals.end(), 0.0);
      }
    }

    // Count nonlinear NZ: nonlinear objective gradient entries + nonlinear constraint entries
    const casadi_int* g_colind_s = jacg_sp_.colind();
    const casadi_int* g_row_s    = jacg_sp_.row();
    int num_nl_nz = 0;
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

    COIDEF_NumCon(m->cntvect, ng_expanded + 1);
    COIDEF_NumNz(m->cntvect, numnz);
    COIDEF_NumNlNz(m->cntvect, num_nl_nz);

    int ret = COI_Solve(m->cntvect);
    if (ret != 0) {
      m->success = false;
      m->unified_return_status = m->nan_encountered ? SOLVER_RET_NAN : SOLVER_RET_UNKNOWN;
      return 0;
    }

    m->success = !m->nan_encountered &&
                 (m->modsta == 1 || m->modsta == 2) && m->solsta == 1;

    if (m->nan_encountered || m->solsta == 5) {
      m->unified_return_status = SOLVER_RET_NAN;
    } else if (m->success) {
      m->unified_return_status = SOLVER_RET_SUCCESS;
    } else {
      if (m->solsta == 2 || m->solsta == 3 ||
          m->solsta == 8 || m->solsta == 15) {
        m->unified_return_status = SOLVER_RET_LIMITED;
      } else if (m->modsta == 4 || m->modsta == 5) {
        m->unified_return_status = SOLVER_RET_INFEASIBLE;
      }
    }
    return 0;
  }

  Dict ConoptInterface::get_stats(void* mem) const {
    Dict stats = Nlpsol::get_stats(mem);
    auto m = static_cast<ConoptMemory*>(mem);
    stats["return_status"] = m->return_status;
    stats["modsta"] = m->modsta;
    stats["solsta"] = m->solsta;
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

    if (opt.second.is_double())     *RVAL = opt.second.to_double();
    else if (opt.second.is_int())   *IVAL = opt.second.to_int();
    else if (opt.second.is_bool())  *LVAL = opt.second.to_bool() ? 1 : 0;

    return 0;
  }

  // --- Progress / Interrupt Callback --- //
  int COI_CALLCONV ConoptInterface::cb_progress(int LEN_INT, const int INTX[], int LEN_RL, const double RL[], const double X[], void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    const ConoptInterface& self = m->self;

    if (!self.fcallback_.is_null()) {
        int phase = (LEN_INT > 1) ? INTX[1] : -1;

        double obj_val = 0.0;
        if (LEN_RL > 1 && phase >= 3) {
            obj_val = RL[1];
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

    // Variable bounds and initial point
    for (int i = 0; i < NUMVAR; ++i) {
      if (!std::isinf(m->d_nlp.lbz[i])) LOWER[i] = m->d_nlp.lbz[i];
      if (!std::isinf(m->d_nlp.ubz[i])) UPPER[i] = m->d_nlp.ubz[i];
      CURR[i] = m->d_nlp.z[i];
    }

    // Constraint types and RHS (row 0 = objective, rows 1..ng_expanded = constraints)
    TYPEX[0] = 3;
    RHS[0]   = 0.0;
    for (int r = 0; r < m->ng_expanded; ++r) {
      TYPEX[r + 1] = m->conopt_type[r];
      RHS[r + 1]   = m->conopt_rhs[r];
    }

    // Jacobian structure — built live from jacg_sp_ with range-row duplication
    const casadi_int* g_colind = self.jacg_sp_.colind();
    const casadi_int* g_row    = self.jacg_sp_.row();
    int nz = 0;
    for (int c = 0; c < NUMVAR; ++c) {
      COLSTA[c] = nz;
      if (self.gradf_col_flag_[c]) {
        casadi_int k = self.gradf_col_to_nz_[c];
        ROWNO[nz]  = 0;
        NLFLAG[nz] = self.gradf_nlflag_[k];
        if (self.gradf_nlflag_[k] == 0) VALUE[nz] = m->gradf_const_vals[k];
        nz++;
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

    try {
        m->arg[0] = m->cached_x.data();
        m->arg[1] = m->d_nlp.p;

        m->res[0] = &m->cached_f;
        m->res[1] = m->cached_grad_f.data();
        self.calc_function(m, "nlp_grad_f");

        m->res[0] = m->cached_g.data();
        m->res[1] = m->cached_jac_g.data();
        self.calc_function(m, "nlp_jac_g");

        m->cache_valid = true;
    } catch (std::exception& ex) {
        casadi::uerr() << ex.what() << std::endl;
        *ERRCNT = 1;
        m->nan_encountered = true;
        m->cache_valid = false;
    } catch (...) {
        *ERRCNT = 1;
        m->nan_encountered = true;
        m->cache_valid = false;
    }

    // Detect silent NaN from CasADi (e.g. sqrt of negative number)
    if (m->cache_valid) {
        bool has_nan = std::isnan(m->cached_f);
        if (!has_nan) {
            for (double v : m->cached_g) if (std::isnan(v)) { has_nan = true; break; }
        }
        if (!has_nan) {
            for (double v : m->cached_grad_f) if (std::isnan(v)) { has_nan = true; break; }
        }
        if (has_nan) {
            *ERRCNT = 1;
            m->nan_encountered = true;
            m->cache_valid = false;
        }
    }
    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_fd_eval(const double X[], double* G, double JAC[], int ROWNO, const int JACNUM[], int MODE, int IGNERR, int* ERRCNT, int NUMVAR, int NUMJAC, int THREAD, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    const ConoptInterface& self = m->self;

    if (!m->cache_valid) {
        *ERRCNT = 1;
        return 0;
    }

    if (ROWNO == 0) {
        if (MODE == 1 || MODE == 3) *G = m->cached_f;
        if (MODE == 2 || MODE == 3) {
            std::memset(JAC, 0, NUMVAR * sizeof(double));
            const casadi_int* f_row = self.gradf_sp_.row();
            for (casadi_int k = 0; k < self.gradf_sp_.nnz(); ++k)
                JAC[f_row[k]] = m->cached_grad_f[k];
        }
    } else {
        int ci = m->conopt_to_casadi[ROWNO - 1];
        if (MODE == 1 || MODE == 3) *G = m->cached_g[ci];
        if (MODE == 2 || MODE == 3) {
            int base = self.jacg_rowstart_[ci];
            for (int k = 0; k < NUMJAC; ++k)
                JAC[JACNUM[k]] = m->cached_jac_g[self.jacg_nzidx_[base + k]];
        }
    }
    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_fdevalend(int IGNERR, int* ERRCNT, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    m->cache_valid = false;
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

    double obj_factor = U[0];
    const double* lam_g = &U[1];

    m->arg[0] = X;
    m->arg[1] = m->d_nlp.p;
    m->arg[2] = &obj_factor;
    m->arg[3] = lam_g;
    m->res[0] = HSVL;

    try {
        self.calc_function(m, "nlp_hess_l");
    } catch (...) {
        *NODRV = 1;
    }
    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_status(int MODSTA, int SOLSTA, int ITER, double OBJVAL, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);
    m->modsta = MODSTA;
    m->solsta = SOLSTA;
    m->iter = ITER;
    m->d_nlp.objective = OBJVAL;

    const char* modsta_str;
    switch (MODSTA) {
      case 1:  modsta_str = "Optimal";                  break;
      case 2:  modsta_str = "Locally optimal";          break;
      case 3:  modsta_str = "Unbounded";                break;
      case 4:  modsta_str = "Infeasible";               break;
      case 5:  modsta_str = "Locally infeasible";       break;
      case 6:  modsta_str = "Intermediate infeasible";  break;
      case 7:  modsta_str = "Intermediate non-optimal"; break;
      case 12: modsta_str = "Unknown error";            break;
      case 13: modsta_str = "Error: no solution";       break;
      default: modsta_str = "Unknown model status";     break;
    }

    const char* solsta_str;
    switch (SOLSTA) {
      case 1:  solsta_str = "Normal completion";                   break;
      case 2:  solsta_str = "Iteration limit";                     break;
      case 3:  solsta_str = "Time limit";                          break;
      case 4:  solsta_str = "Terminated by solver";                break;
      case 5:  solsta_str = "Evaluation error limit";              break;
      case 8:  solsta_str = "User interrupt";                      break;
      case 9:  solsta_str = "Setup failure";                       break;
      case 10: solsta_str = "Major solver error";                  break;
      case 11: solsta_str = "Major solver error (feasible point)"; break;
      case 13: solsta_str = "System error";                        break;
      case 15: solsta_str = "Quick Mode termination";              break;
      default: solsta_str = "Unknown solver status";               break;
    }

    m->return_status = std::string(modsta_str) + " / " + std::string(solsta_str);
    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_solution(const double XVAL[], const double XMAR[], const int XBAS[], const int XSTA[], const double YVAL[], const double YMAR[], const int YBAS[], const int YSTA[], int NUMVAR, int NUMCON, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);

    casadi_copy(XVAL, NUMVAR, m->d_nlp.z);

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
    for (int i = 0; i < SMSG; ++i) {
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
