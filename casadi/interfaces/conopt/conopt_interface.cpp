#include "conopt_interface.hpp"
#include "casadi/core/casadi_misc.hpp"
#include "casadi/core/casadi_interrupt.hpp"
#include <cmath>
#include <cstring>
#include <algorithm>

namespace casadi {

  extern "C" int casadi_register_nlpsol_conopt(Nlpsol::Plugin* plugin) {
    plugin->creator = ConoptInterface::creator;
    plugin->name = "conopt";
    plugin->doc = "CONOPT Interface";
    plugin->version = CASADI_VERSION;
    plugin->options = &ConoptInterface::options_;
    plugin->deserialize = &ConoptInterface::deserialize;
    return 0;
  }

  extern "C" void casadi_load_nlpsol_conopt() {
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

    // Setup 2nd Order Info
    exact_hessian_ = true;
    if (opts.find("exact_hessian") != opts.end()) exact_hessian_ = opts.at("exact_hessian");

    if (exact_hessian_) {
      Function hl_fcn = create_function("nlp_hess_l", {"x", "p", "lam:f", "lam:g"},
                                        {"tril:hess:gamma:x:x"}, {{"gamma", {"f", "g"}}});
      hesslag_sp_ = hl_fcn.sparsity_out(0);
      alloc_w(hesslag_sp_.nnz(), true);
    }

    jac_colsta_.resize(nx_ + 1, 0);
    int numnz = 0;
    std::vector<bool> has_gradf(nx_, false);
    const casadi_int* f_row = gradf_sp_.row();
    for (casadi_int k = 0; k < gradf_sp_.nnz(); ++k) has_gradf[f_row[k]] = true;

    const casadi_int* g_colind = jacg_sp_.colind();
    const casadi_int* g_row = jacg_sp_.row();

    for (int c = 0; c < nx_; ++c) {
      jac_colsta_[c] = numnz;
      if (has_gradf[c]) {
        jac_rowno_.push_back(0);
        jac_nlflag_.push_back(1);
        numnz++;
      }
      for (casadi_int el = g_colind[c]; el < g_colind[c+1]; ++el) {
        jac_rowno_.push_back(g_row[el] + 1);
        jac_nlflag_.push_back(1);
        numnz++;
      }
    }
    jac_colsta_[nx_] = numnz;

    alloc_w(1, true);
    alloc_w(gradf_sp_.nnz(), true);
    alloc_w(ng_, true);
    alloc_w(jacg_sp_.nnz(), true);
  }

  // --- Serialization & Deserialization --- //
  ConoptInterface::ConoptInterface(DeserializingStream& s) : Nlpsol(s) {
    int version = s.version("ConoptInterface", 1);
    s.unpack("ConoptInterface::exact_hessian", exact_hessian_);
    s.unpack("ConoptInterface::opts", opts_);
    s.unpack("ConoptInterface::gradf_sp", gradf_sp_);
    s.unpack("ConoptInterface::jacg_sp", jacg_sp_);
    s.unpack("ConoptInterface::hesslag_sp", hesslag_sp_);
    s.unpack("ConoptInterface::jac_colsta", jac_colsta_);
    s.unpack("ConoptInterface::jac_rowno", jac_rowno_);
    s.unpack("ConoptInterface::jac_nlflag", jac_nlflag_);
  }

  void ConoptInterface::serialize_body(SerializingStream &s) const {
    Nlpsol::serialize_body(s);
    s.version("ConoptInterface", 1);
    s.pack("ConoptInterface::exact_hessian", exact_hessian_);
    s.pack("ConoptInterface::opts", opts_);
    s.pack("ConoptInterface::gradf_sp", gradf_sp_);
    s.pack("ConoptInterface::jacg_sp", jacg_sp_);
    s.pack("ConoptInterface::hesslag_sp", hesslag_sp_);
    s.pack("ConoptInterface::jac_colsta", jac_colsta_);
    s.pack("ConoptInterface::jac_rowno", jac_rowno_);
    s.pack("ConoptInterface::jac_nlflag", jac_nlflag_);
  }

  ConoptMemory::ConoptMemory(const ConoptInterface& interface)
      : self(interface), NlpsolMemory(), cntvect(nullptr), modsta(0), solsta(0),
        iter(0), return_status("Unset"), success(0), cache_valid(false), current_option_idx(0) {}

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

    COI_Create(&m->cntvect);

    COIDEF_NumVar(m->cntvect, nx_);
    COIDEF_NumCon(m->cntvect, ng_ + 1);
    COIDEF_NumNz(m->cntvect, jac_rowno_.size());

    COIDEF_ObjCon(m->cntvect, 0);
    COIDEF_OptDir(m->cntvect, 1);

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
    m->current_option_idx = 0;
    COIDEF_Option(m->cntvect, &ConoptInterface::cb_option);
    COIDEF_Progress(m->cntvect, &ConoptInterface::cb_progress);

    // Register Callbacks
    COIDEF_ReadMatrix(m->cntvect, &ConoptInterface::cb_read_matrix);
    COIDEF_FDEvalIni(m->cntvect, &ConoptInterface::cb_fdevalini);
    COIDEF_FDEval(m->cntvect, &ConoptInterface::cb_fd_eval);
    COIDEF_FDEvalEnd(m->cntvect, &ConoptInterface::cb_fdevalend);

    if (exact_hessian_) {
        COIDEF_NumHess(m->cntvect, hesslag_sp_.nnz());
        COIDEF_2DLagrSize(m->cntvect, &ConoptInterface::cb_2dlagrsize);
        COIDEF_2DLagrStr(m->cntvect, &ConoptInterface::cb_2dlagrstr);
        COIDEF_2DLagrVal(m->cntvect, &ConoptInterface::cb_2dlagrval);
    }

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
    COI_Solve(m->cntvect);
    m->success = (m->solsta == 1);
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

    if (m->current_option_idx < m->custom_options.size()) {
        auto& opt = m->custom_options[m->current_option_idx];

        std::string name = opt.first;
        if (name.length() > 8) name = name.substr(0, 8);
        else name.append(8 - name.length(), ' ');
        std::strncpy(NAME, name.c_str(), 8);

        if (opt.second.is_double()) *RVAL = opt.second.to_double();
        else if (opt.second.is_int()) *IVAL = opt.second.to_int();
        else if (opt.second.is_bool()) *LVAL = opt.second.to_bool() ? 1 : 0;

        m->current_option_idx++;
    } else {
        std::strncpy(NAME, "        ", 8);
    }
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

    for(int i = 0; i < NUMVAR; ++i) {
      if (!std::isinf(m->d_nlp.lbz[i])) LOWER[i] = m->d_nlp.lbz[i];
      if (!std::isinf(m->d_nlp.ubz[i])) UPPER[i] = m->d_nlp.ubz[i];
      CURR[i]  = m->d_nlp.z[i];
    }

    TYPEX[0] = 3;
    RHS[0] = 0.0;

    for(int i = 0; i < NUMCON - 1; ++i) {
      double lbg = m->d_nlp.lbz[NUMVAR + i];
      double ubg = m->d_nlp.ubz[NUMVAR + i];
      int row_idx = i + 1;

      if (lbg == ubg) {
        TYPEX[row_idx] = 0;
        RHS[row_idx] = lbg;
      } else if (!std::isinf(lbg) && std::isinf(ubg)) {
        TYPEX[row_idx] = 1;
        RHS[row_idx] = lbg;
      } else if (std::isinf(lbg) && !std::isinf(ubg)) {
        TYPEX[row_idx] = 2;
        RHS[row_idx] = ubg;
      } else if (std::isinf(lbg) && std::isinf(ubg)) {
        TYPEX[row_idx] = 3;
        RHS[row_idx] = 0.0;
      } else {
        casadi_error("CONOPT does not accept range constraints.");
      }
    }

    std::memcpy(COLSTA, self.jac_colsta_.data(), (NUMVAR + 1) * sizeof(int));
    std::memcpy(ROWNO, self.jac_rowno_.data(), NUMNZ * sizeof(int));
    std::memcpy(NLFLAG, self.jac_nlflag_.data(), NUMNZ * sizeof(int));

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
    } catch (...) {
        *ERRCNT = 1;
        m->cache_valid = false;
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
            for (casadi_int k = 0; k < self.gradf_sp_.nnz(); ++k) {
                JAC[f_row[k]] = m->cached_grad_f[k];
            }
        }
    }
    else {
        int g_rowno = ROWNO - 1;
        if (MODE == 1 || MODE == 3) *G = m->cached_g[g_rowno];
        if (MODE == 2 || MODE == 3) {
            std::memset(JAC, 0, NUMVAR * sizeof(double));
            const casadi_int* colind = self.jacg_sp_.colind();
            const casadi_int* row = self.jacg_sp_.row();
            for (int c = 0; c < NUMVAR; ++c) {
                for (casadi_int el = colind[c]; el < colind[c+1]; ++el) {
                    if (row[el] == g_rowno) {
                        JAC[c] = m->cached_jac_g[el];
                    }
                }
            }
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

    if (SOLSTA == 1) m->return_status = "Optimal";
    else if (SOLSTA >= 6 && SOLSTA <= 11) m->return_status = "No Solution (Major Error)";
    else m->return_status = "Suboptimal/Infeasible";

    return 0;
  }

  int COI_CALLCONV ConoptInterface::cb_solution(const double XVAL[], const double XMAR[], const int XBAS[], const int XSTA[], const double YVAL[], const double YMAR[], const int YBAS[], const int YSTA[], int NUMVAR, int NUMCON, void* USRMEM) {
    auto m = static_cast<ConoptMemory*>(USRMEM);

    casadi_copy(XVAL, NUMVAR, m->d_nlp.z);
    casadi_copy(XMAR, NUMVAR, m->d_nlp.lam);
    casadi_copy(YVAL + 1, NUMCON - 1, m->d_nlp.z + NUMVAR);
    casadi_copy(YMAR + 1, NUMCON - 1, m->d_nlp.lam + NUMVAR);

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
