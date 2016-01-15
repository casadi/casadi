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


#include "worhp_interface.hpp"

#include "casadi/core/std_vector_tools.hpp"
#include <ctime>
#include <cstring>

using namespace std;

namespace casadi {

  extern "C"
  int CASADI_NLPSOL_WORHP_EXPORT
  casadi_register_nlpsol_worhp(Nlpsol::Plugin* plugin) {
    plugin->creator = WorhpInterface::creator;
    plugin->name = "worhp";
    plugin->doc = WorhpInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_WORHP_EXPORT casadi_load_nlpsol_worhp() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_worhp);
  }

  WorhpInterface::WorhpInterface(const std::string& name, const XProblem& nlp)
    : Nlpsol(name, nlp) {
  }

  WorhpInterface::~WorhpInterface() {
  }

  Options WorhpInterface::options_
  = {{&Nlpsol::options_},
     {{"worhp",
       {OT_DICT,
        "Options to be passed to WORHP"}},
      {"print_time",
       {OT_BOOL,
        "Print information about execution time"}}
     }
  };

  void WorhpInterface::init(const Dict& opts) {

    // Call the init method of the base class
    Nlpsol::init(opts);

    // Default options
    Dict worhp_opts;
    print_time_ = true;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="worhp") {
        worhp_opts = op.second;
      } else if (op.first=="print_time") {
        print_time_ = op.second;
      }
    }

    // Sort Worhp options
    int nopts = WorhpGetParamCount();
    for (auto&& op : worhp_opts) {
      // Get corresponding index using a linear search
      int ind;
      for (ind=1; ind<=nopts; ++ind) {
        // Get name in WORHP
        const char* name = WorhpGetParamName(ind);
        // Break if matching name
        if (op.first.compare(name)==0) break;
      }
      if (ind>nopts) casadi_error("No such Worhp option: " + op.first);

      // Add to the corresponding list
      switch (WorhpGetParamType(ind)) {
      case WORHP_BOOL_T:
        bool_opts_[op.first] = op.second;
        break;
      case WORHP_DOUBLE_T:
        double_opts_[op.first] = op.second;
        break;
      case WORHP_INT_T:
        int_opts_[op.first] = op.second;
        break;
      default:
        casadi_error("Cannot handle WORHP option \"" + op.first + "\": Unknown type");
        break;
      }
    }

    // Setup NLP functions
    setup_f(); // Objective
    setup_g(); // Constraints
    setup_grad_f(); // Gradient of objective
    setup_jac_g(); // Jacobian of the constraints
    setup_hess_l(true); // Hessian of the Lagrangian

    // Temporary vectors
    alloc_w(nx_); // for fetching diagonal entries form Hessian
  }

  void WorhpInterface::init_memory(Memory& mem) const {
    Nlpsol::init_memory(mem);
    WorhpMemory& m = dynamic_cast<WorhpMemory&>(mem);

    // Initialize parameters to default values
    int status;
    InitParams(&status, &m.worhp_p);

    // Pass boolean parameters
    for (auto&& op : bool_opts_) {
      WorhpSetBoolParam(&m.worhp_p, op.first.c_str(), op.second);
    }

    // Pass double parameters
    for (auto&& op : double_opts_) {
      WorhpSetDoubleParam(&m.worhp_p, op.first.c_str(), op.second);
    }

    // Pass integer parameters
    for (auto&& op : int_opts_) {
      WorhpSetIntParam(&m.worhp_p, op.first.c_str(), op.second);
    }

    // Mark the parameters as set
    m.worhp_p.initialised = true;
  }

  void WorhpInterface::set_work(Memory& mem, const double**& arg, double**& res,
                                int*& iw, double*& w) const {
    WorhpMemory& m = dynamic_cast<WorhpMemory&>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    // Number of (free) variables
    m.worhp_o.n = nx_;

    // Number of constraints
    m.worhp_o.m = ng_;

    // Free existing Worhp memory (except parameters)
    bool p_init_backup = m.worhp_p.initialised;
    m.worhp_p.initialised = false; // Avoid freeing the memory for parameters
    if (m.worhp_o.initialised || m.worhp_w.initialised || m.worhp_c.initialised) {
      WorhpFree(&m.worhp_o, &m.worhp_w, &m.worhp_p, &m.worhp_c);
    }
    m.worhp_p.initialised = p_init_backup;

    /// Control data structure needs to be reset every time
    m.worhp_c.initialised = false;
    m.worhp_w.initialised = false;
    m.worhp_o.initialised = false;

    // Worhp uses the CS format internally, hence it is the preferred sparse matrix format.
    m.worhp_w.DF.nnz = nx_;
    if (m.worhp_o.m>0) {
      m.worhp_w.DG.nnz = jacg_sp_.nnz();  // Jacobian of G
    } else {
      m.worhp_w.DG.nnz = 0;
    }

    if (true /*m.worhp_w.HM.NeedStructure*/) { // not initialized
      m.worhp_w.HM.nnz = nx_ + hesslag_sp_.nnz_lower(true);
    } else {
      m.worhp_w.HM.nnz = 0;
    }

    /* Data structure initialisation. */
    WorhpInit(&m.worhp_o, &m.worhp_w, &m.worhp_p, &m.worhp_c);
    if (m.worhp_c.status != FirstCall) {
      string msg = return_codes(m.worhp_c.status);
      casadi_error("Main: Initialisation failed. Status: " + msg);
    }

    if (m.worhp_w.DF.NeedStructure) {
      for (int i=0; i<nx_; ++i) {
        m.worhp_w.DF.row[i] = i + 1; // Index-1 based
      }
    }

    if (m.worhp_o.m>0 && m.worhp_w.DG.NeedStructure) {
      int nz=0;
      const int* colind = jacg_sp_.colind();
      const int* row = jacg_sp_.row();
      for (int c=0; c<nx_; ++c) {
        for (int el=colind[c]; el<colind[c+1]; ++el) {
          int r = row[el];
          m.worhp_w.DG.col[nz] = c + 1; // Index-1 based
          m.worhp_w.DG.row[nz] = r + 1;
          nz++;
        }
      }
    }

    if (m.worhp_w.HM.NeedStructure) {
      // Get the sparsity pattern of the Hessian
      const int* colind = hesslag_sp_.colind();
      const int* row = hesslag_sp_.row();

      int nz=0;

      // Strictly lower triangular part of the Hessian (note CCS -> CRS format change)
      for (int c=0; c<nx_; ++c) {
        for (int el=colind[c]; el<colind[c+1]; ++el) {
          if (row[el]>c) {
            m.worhp_w.HM.row[nz] = row[el] + 1;
            m.worhp_w.HM.col[nz] = c + 1;
            nz++;
          }
        }
      }

      // Diagonal always included
      for (int r=0; r<nx_; ++r) {
        m.worhp_w.HM.row[nz] = r + 1;
        m.worhp_w.HM.col[nz] = r + 1;
        nz++;
      }
    }
  }

  void WorhpInterface::solve(Memory& mem) const {
    WorhpMemory& m = dynamic_cast<WorhpMemory&>(mem);

    // Check the provided inputs
    checkInputs(mem);

    // Reset the counters
    m.t_eval_f = m.t_eval_grad_f = m.t_eval_g = m.t_eval_jac_g = m.t_eval_h = m.t_callback_fun =
      m.t_callback_prepare = m.t_mainloop = 0;
    m.n_eval_f = m.n_eval_grad_f = m.n_eval_g = m.n_eval_jac_g = m.n_eval_h = 0;

    double inf = numeric_limits<double>::infinity();

    if (m.lbx && m.ubx) {
      for (int i=0; i<nx_;++i) {
        casadi_assert_message(m.lbx[i]!=m.ubx[i],
                              "WorhpInterface::evaluate: Worhp cannot handle the case when "
                              "LBX == UBX."
                              "You have that case at non-zero " << i << " , which has value " <<
                              m.ubx[i] << ". Reformulate your problem by using a parameter "
                              "for the corresponding variable.");
      }
    }

    if (m.lbg && m.ubg) {
      for (int i=0; i<ng_; ++i) {
        casadi_assert_message(!(m.lbg[i]==-inf && m.ubg[i] == inf),
                              "WorhpInterface::evaluate: Worhp cannot handle the case when both "
                              "LBG and UBG are infinite."
                              "You have that case at non-zero " << i << "."
                              "Reformulate your problem eliminating the corresponding constraint.");
      }
    }

    // Pass inputs to WORHP data structures
    casadi_copy(m.x0, nx_, m.worhp_o.X);
    casadi_copy(m.lbx, nx_, m.worhp_o.XL);
    casadi_copy(m.ubx, nx_, m.worhp_o.XU);
    casadi_copy(m.lam_x0, nx_, m.worhp_o.Lambda);
    if (m.worhp_o.m>0) {
      casadi_copy(m.lam_g0, ng_, m.worhp_o.Mu);
      casadi_copy(m.lbg, ng_, m.worhp_o.GL);
      casadi_copy(m.ubg, ng_, m.worhp_o.GU);
    }

    // Replace infinite bounds with m.worhp_p.Infty
    for (int i=0; i<nx_; ++i) if (m.worhp_o.XL[i]==-inf) m.worhp_o.XL[i] = -m.worhp_p.Infty;
    for (int i=0; i<nx_; ++i) if (m.worhp_o.XU[i]== inf) m.worhp_o.XU[i] =  m.worhp_p.Infty;
    for (int i=0; i<ng_; ++i) if (m.worhp_o.GL[i]==-inf) m.worhp_o.GL[i] = -m.worhp_p.Infty;
    for (int i=0; i<ng_; ++i) if (m.worhp_o.GU[i]== inf) m.worhp_o.GU[i] =  m.worhp_p.Infty;

    log("WorhpInterface::starting iteration");

    double time1 = clock();

    bool firstIteration = true;

    // Reverse Communication loop
    while (m.worhp_c.status < TerminateSuccess &&  m.worhp_c.status > TerminateError) {
      if (GetUserAction(&m.worhp_c, callWorhp)) {
        Worhp(&m.worhp_o, &m.worhp_w, &m.worhp_p, &m.worhp_c);
      }


      if (GetUserAction(&m.worhp_c, iterOutput)) {

        if (!firstIteration) {
          firstIteration = true;

          if (!fcallback_.isNull()) {
            m.iter = m.worhp_w.MajorIter;
            m.iter_sqp = m.worhp_w.MinorIter;
            m.inf_pr = m.worhp_w.NormMax_CV;
            m.inf_du = m.worhp_w.ScaledKKT;
            m.alpha_pr = m.worhp_w.ArmijoAlpha;

            time1 = clock();

            // Inputs
            fill_n(m.arg, fcallback_.n_in(), nullptr);
            m.arg[NLPSOL_X] = m.worhp_o.X;
            m.arg[NLPSOL_F] = &m.worhp_o.F;
            m.arg[NLPSOL_G] = m.worhp_o.G;
            m.arg[NLPSOL_LAM_P] = 0;
            m.arg[NLPSOL_LAM_X] = m.worhp_o.Lambda;
            m.arg[NLPSOL_LAM_G] = m.worhp_o.Mu;

            // Outputs
            fill_n(m.res, fcallback_.n_out(), nullptr);
            double ret_double;
            m.res[0] = &ret_double;

            // Evaluate the callback function
            m.n_eval_callback += 1;
            fcallback_(m.arg, m.res, m.iw, m.w, 0);
            int ret = static_cast<int>(ret_double);

            if (ret) m.worhp_c.status = TerminatedByUser;
          }
        }


        IterationOutput(&m.worhp_o, &m.worhp_w, &m.worhp_p, &m.worhp_c);
        DoneUserAction(&m.worhp_c, iterOutput);
      }

      if (GetUserAction(&m.worhp_c, evalF)) {
        calc_f(m, m.worhp_o.X, m.p, &m.worhp_o.F);
        if (m.f) *m.f = m.worhp_o.F; // Store cost, before scaling
        m.worhp_o.F *= m.worhp_w.ScaleObj;
        DoneUserAction(&m.worhp_c, evalF);
      }

      if (GetUserAction(&m.worhp_c, evalG)) {
        calc_g(m, m.worhp_o.X, m.p, m.worhp_o.G);
        DoneUserAction(&m.worhp_c, evalG);
      }

      if (GetUserAction(&m.worhp_c, evalDF)) {
        calc_grad_f(m, m.worhp_o.X, m.p, 0, m.worhp_w.DF.val);
        casadi_scal(nx_, m.worhp_w.ScaleObj, m.worhp_w.DF.val);
        DoneUserAction(&m.worhp_c, evalDF);
      }

      if (GetUserAction(&m.worhp_c, evalDG)) {
        calc_jac_g(m, m.worhp_o.X, m.p, 0, m.worhp_w.DG.val);
        DoneUserAction(&m.worhp_c, evalDG);
      }

      if (GetUserAction(&m.worhp_c, evalHM)) {
        calc_hess_l(m, m.worhp_o.X, m.p, &m.worhp_w.ScaleObj, m.worhp_o.Mu,
                    m.worhp_w.HM.val);
        // Diagonal values
        double *dval = m.w;
        casadi_fill(dval, nx_, 0.);

        // Remove diagonal
        const int* colind = hesslag_sp_.colind();
        const int* row = hesslag_sp_.row();
        int ind=0;
        for (int c=0; c<nx_; ++c) {
          for (int el=colind[c]; el<colind[c+1]; ++el) {
            if (row[el]==c) {
              dval[c] = m.worhp_w.HM.val[el];
            } else {
              m.worhp_w.HM.val[ind++] = m.worhp_w.HM.val[el];
            }
          }
        }

        // Add diagonal entries at the end
        casadi_copy(dval, nx_, m.worhp_w.HM.val+ind);
        DoneUserAction(&m.worhp_c, evalHM);
      }

      if (GetUserAction(&m.worhp_c, fidif)) {
        WorhpFidif(&m.worhp_o, &m.worhp_w, &m.worhp_p, &m.worhp_c);
      }
    }

    double time2 = clock();
    m.t_mainloop += (time2-time1)/CLOCKS_PER_SEC;

    // Copy outputs
    casadi_copy(m.worhp_o.X, nx_, m.x);
    casadi_copy(m.worhp_o.G, ng_, m.g);
    casadi_copy(m.worhp_o.Lambda, nx_, m.lam_x);
    casadi_copy(m.worhp_o.Mu, ng_, m.lam_g);

    StatusMsg(&m.worhp_o, &m.worhp_w, &m.worhp_p, &m.worhp_c);

    if (print_time_) {
      // Write timings
      userOut() << "time spent in eval_f: " << m.t_eval_f << " s.";
      if (m.n_eval_f>0)
        userOut() << " (" << m.n_eval_f << " calls, " <<
          (m.t_eval_f/m.n_eval_f)*1000 << " ms. average)";
      userOut() << endl;
      userOut() << "time spent in eval_grad_f: " << m.t_eval_grad_f << " s.";
      if (m.n_eval_grad_f>0)
        userOut() << " (" << m.n_eval_grad_f << " calls, "
             << (m.t_eval_grad_f/m.n_eval_grad_f)*1000 << " ms. average)";
      userOut() << endl;
      userOut() << "time spent in eval_g: " << m.t_eval_g << " s.";
      if (m.n_eval_g>0)
        userOut() << " (" << m.n_eval_g << " calls, " <<
          (m.t_eval_g/m.n_eval_g)*1000 << " ms. average)";
      userOut() << endl;
      userOut() << "time spent in eval_jac_g: " << m.t_eval_jac_g << " s.";
      if (m.n_eval_jac_g>0)
        userOut() << " (" << m.n_eval_jac_g << " calls, "
                  << (m.t_eval_jac_g/m.n_eval_jac_g)*1000 << " ms. average)";
      userOut() << endl;
      userOut() << "time spent in eval_h: " << m.t_eval_h << " s.";
      if (m.n_eval_h>1)
        userOut() << " (" << m.n_eval_h << " calls, " <<
          (m.t_eval_h/m.n_eval_h)*1000 << " ms. average)";
      userOut() << endl;
      userOut() << "time spent in main loop: " << m.t_mainloop << " s." << endl;
      userOut() << "time spent in callback function: " << m.t_callback_fun << " s." << endl;
      userOut() << "time spent in callback preparation: " << m.t_callback_prepare << " s." << endl;
    }

    m.return_code = m.worhp_c.status;
    m.return_status = return_codes(m.worhp_c.status);
  }

  const char* WorhpInterface::return_codes(int flag) {
    switch (flag) {
    case TerminateSuccess: return "TerminateSuccess";
    case OptimalSolution: return "OptimalSolution";
    case SearchDirectionZero: return "SearchDirectionZero";
    case SearchDirectionSmall: return "SearchDirectionSmall";
    case StationaryPointFound: return "StationaryPointFound";
    case AcceptablePrevious: return "AcceptablePrevious";
    case FritzJohn: return "FritzJohn";
    case NotDiffable: return "NotDiffable";
    case Unbounded: return "Unbounded";
    case FeasibleSolution: return "FeasibleSolution";
    case LowPassFilterOptimal: return "LowPassFilterOptimal";
    case LowPassFilterAcceptable: return "LowPassFilterAcceptable";
    case TerminateError: return "TerminateError";
    case InitError: return "InitError";
    case DataError: return "DataError";
    case MaxCalls: return "MaxCalls";
    case MaxIter: return "MaxIter";
    case MinimumStepsize: return "MinimumStepsize";
    case QPerror: return "QPerror";
    case ProblemInfeasible: return "ProblemInfeasible";
    case GroupsComposition: return "GroupsComposition";
    case TooBig: return "TooBig";
    case Timeout: return "Timeout";
    case FDError: return "FDError";
    case LocalInfeas: return "LocalInfeas";
    case LicenseError: return "LicenseError";
    case TerminatedByUser: return "TerminatedByUser";
    case FunctionErrorF: return "FunctionErrorF";
    case FunctionErrorG: return "FunctionErrorG";
    case FunctionErrorDF: return "FunctionErrorDF";
    case FunctionErrorDG: return "FunctionErrorDG";
    case FunctionErrorHM: return "FunctionErrorHM";
    }
    return "Unknown WORHP return code";
  }

  WorhpMemory::WorhpMemory() {
    this->worhp_o.initialised = false;
    this->worhp_w.initialised = false;
    this->worhp_p.initialised = false;
    this->worhp_c.initialised = false;
  }

  WorhpMemory::~WorhpMemory() {
    if (this->worhp_p.initialised || this->worhp_o.initialised ||
        this->worhp_w.initialised || this->worhp_c.initialised) {
      WorhpFree(&this->worhp_o, &this->worhp_w, &this->worhp_p, &this->worhp_c);
    }
  }

} // namespace casadi
