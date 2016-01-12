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


#include "knitro_interface.hpp"
#include "casadi/core/std_vector_tools.hpp"
#include <ctime>
#include <cstdio>
#include <cstdlib>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_NLPSOL_KNITRO_EXPORT
  casadi_register_nlpsol_knitro(Nlpsol::Plugin* plugin) {
    plugin->creator = KnitroInterface::creator;
    plugin->name = "knitro";
    plugin->doc = KnitroInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_KNITRO_EXPORT casadi_load_nlpsol_knitro() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_knitro);
  }

  KnitroInterface::KnitroInterface(const std::string& name, const XProblem& nlp)
    : Nlpsol(name, nlp) {

    // Monitors
    addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "",
              "eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h", true);

    // Not yet ready
    //addOption("algorithm",                OT_STRING, GenericType(),
    // "Which algorithm to use. See KNITRO documentation.", "auto|direct|cg|active");
    //addOption("bar_directinterval",       OT_INT, GenericType(),
    //  "When using the Interior/Direct algorithm, this parameter controls the maximum number of "
    //  "consecutive CG steps before trying to force the algorithm to take a direct step again. "
    //  "See KNITRO documentation.");
    //addOption("bar_feasible",             OT_STRING, GenericType(),
    //  "Whether feasibility is given special emphasis. See KNITRO documentation.",
    //  "no|stay|get|get_stay");
    //addOption("bar_feasmodetol",          OT_DOUBLE, GenericType(),
    //  "Specifies the tolerance for entering the stay feasible mode See KNITRO documentation.");
    //addOption("bar_initmu",               OT_INT, GenericType(),
    //  "Initial value for the barrier parameter. See KNITRO documentation.");
    //addOption("bar_initpt",               OT_STRING, GenericType(),
    //  "Whether to use the initial point strategy with barrier algorithms. "
    //  "See KNITRO documentation.", "auto|yes|no");
    //addOption("bar_maxbacktrack",         OT_INT, GenericType(),
    //  "Maximum allowable number of backtracks during the linesearch of the Interior Direct "
    //  "algorithm before reverting to a CG step. See KNITRO documentation.");
    //addOption("bar_maxrefactor",          OT_INT, GenericType(),
    //  "Maximum number of refactorizations of the KKT system per iteration of the Interior "
    //  "Direct algorithm before reverting to a CG step. See KNITRO documentation.");

    //addOption("Alg", OT_INT,0, "Algorithm");
    addOption("BarRule", OT_INT, 0, "Barrier Rule");
    addOption("NewPoint", OT_BOOL, 0, "Select new-point feature");
    addOption("GradOpt", OT_INT, 1, "Gradient calculation method");
    addOption("HessOpt", OT_INT, 1, "Hessian calculation method");
    addOption("Feasible", OT_BOOL, 1, "Allow infeasible iterations");
    addOption("HonorBnds", OT_BOOL, 0, "Enforce bounds");
    addOption("LpSolver", OT_BOOL, 0, "Use LpSolver");
    addOption("Multistart", OT_BOOL, 0, "Use multistart");
    //addOption("MsMaxSolves", OT_INT, 1, "Maximum multistart points");
    addOption("MaxCgIt", OT_INT, 0, "Maximum conjugate gradient iterations");
    //addOption("MaxCrossTt", OT_INT, 0, "Maximum crossover iterations");
    addOption("MaxIt", OT_INT, 10000, "Iteration limit");
    //addOption("MaxTimeCPU", OT_DOUBLE, 1e8, "CPU Time limit");
    //addOption("MaxTimeReal", OT_DOUBLE, 1e8, "Time limit");
    addOption("LmSize", OT_INT, 10, "Memory pairsize limit");
    addOption("Scale", OT_BOOL, 1, "Perform scaling");
    addOption("ShiftInit", OT_BOOL, 1, "Interior-point shifting initial point");
    addOption("Soc", OT_INT, 1, "Second order correction");
    addOption("InitPt", OT_BOOL, 0, "Use initial point strategy");
    addOption("Delta", OT_DOUBLE, 1.0, "Initial region scaling factor");
    addOption("FeasModeTol", OT_DOUBLE, 1e-4, "Feasible mode tolerance");
    addOption("FeasTol", OT_DOUBLE, 1e-6, "Feasible tolerance");
    addOption("FeasTolAbs", OT_DOUBLE, 0, "Absolute feasible tolerance");
    addOption("OptTol", OT_DOUBLE, 1e-6, "Relative optimality tolerance");
    addOption("OptTolAbs", OT_DOUBLE, 0, "Absolute optimality tolerance");
    addOption("Pivot", OT_DOUBLE, 1e-8, "Initial pivot threshold");
    addOption("XTol", OT_DOUBLE, 1e-15, "Relative solution change tolerance");
    addOption("Mu", OT_DOUBLE, 0.1, "Initial barrier parameter");
    addOption("ObjRange", OT_DOUBLE, 1e20, "Maximum objective value");
    addOption("OutLev", OT_INT, 2, "Log output level");
    addOption("Debug", OT_INT, 0, "Debug level");
    addOption("contype", OT_INTVECTOR);
  }


  KnitroInterface::~KnitroInterface() {
  }

  void KnitroInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    //if (hasSetOption("Alg")) int_param_["alg"] = option("Alg");
    if (hasSetOption("BarRule")) int_param_["barrule"] = option("BarRule");
    if (hasSetOption("NewPoint")) int_param_["newpoint"] = option("NewPoint");
    if (hasSetOption("GradOpt")) int_param_["gradopt"] = option("GradOpt");
    if (hasSetOption("HessOpt")) int_param_["hessopt"] = option("HessOpt");
    if (hasSetOption("Feasible")) int_param_["feasible"] = option("Feasible");
    if (hasSetOption("HonorBnds")) int_param_["honorbnds"] = option("HonorBnds");
    if (hasSetOption("LpSolver")) int_param_["lpsolver"] = option("LpSolver");
    if (hasSetOption("Multistart")) int_param_["multistart"] = option("Multistart");
    //if (hasSetOption("MsMaxSolves")) int_param_["msmaxsolves"] = option("MsMaxSolves");
    if (hasSetOption("MaxCgIt")) int_param_["maxcgit"] = option("MaxCgIt");
    //if (hasSetOption("MaxCrossTt")) int_param_["maxcrosstt"] = option("MaxCrossTt");
    if (hasSetOption("MaxIt")) int_param_["maxit"] = option("MaxIt");
    //if (hasSetOption("MaxTimeCPU")) double_param_["maxtimecpu"] = option("MaxTimeCPU");
    //if (hasSetOption("MaxTimeReal")) double_param_["maxtimereal"] = option("MaxTimeReal");
    if (hasSetOption("LmSize")) int_param_["lmsize"] = option("LmSize");
    if (hasSetOption("Scale")) int_param_["scale"] = option("Scale");
    if (hasSetOption("ShiftInit")) int_param_["shiftinit"] = option("ShiftInit");
    if (hasSetOption("Soc")) int_param_["soc"] = option("Soc");
    if (hasSetOption("InitPt")) int_param_["initpt"] = option("InitPt");
    if (hasSetOption("Delta")) double_param_["delta"] = option("Delta");
    if (hasSetOption("FeasModeTol")) double_param_["feasmodetol"] = option("FeasModeTol");
    if (hasSetOption("FeasTol")) double_param_["feastol"] = option("FeasTol");
    if (hasSetOption("FeasTolAbs")) double_param_["feastolabs"] = option("FeasTolAbs");
    if (hasSetOption("OptTol")) double_param_["opttol"] = option("OptTol");
    if (hasSetOption("OptTolAbs")) double_param_["opttolabs"] = option("OptTolAbs");
    if (hasSetOption("Pivot")) double_param_["pivot"] = option("Pivot");
    if (hasSetOption("XTol")) double_param_["xtol"] = option("XTol");
    if (hasSetOption("Mu")) double_param_["mu"] = option("Mu");
    if (hasSetOption("ObjRange")) double_param_["objrange"] = option("ObjRange");
    if (hasSetOption("OutLev")) int_param_["outlev"] = option("OutLev");
    if (hasSetOption("Debug")) int_param_["debug"] = option("Debug");

    // Setup NLP functions
    setup_fg(); // Objective and constraints
    setup_gf_jg(); // Objective gradient and Jacobian of constraints
    setup_hess_l(); // Hessian of the Lagrangian

    // Allocate memory
    alloc_w(nx_, true); // wx_
    alloc_w(nx_, true); // wlbx_
    alloc_w(nx_, true); // wubx_
    alloc_w(ng_, true); // wlbg_
    alloc_w(ng_, true); // wubg_
  }

  void KnitroInterface::init_memory(Memory& mem) const {
    Nlpsol::init_memory(mem);
    KnitroMemory& m = dynamic_cast<KnitroMemory&>(mem);

    // Commented out since I have not found out how to change the bounds
    // Allocate KNITRO memory block
    /*  m.kc_handle = KTR_new(); */
  }

  void KnitroInterface::set_work(Memory& mem, const double**& arg, double**& res,
                                 int*& iw, double*& w) const {
    KnitroMemory& m = dynamic_cast<KnitroMemory&>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    // Copy inputs to temporary arrays
    m.wx = w; w += nx_;
    m.wlbx = w; w += nx_;
    m.wubx = w; w += nx_;
    m.wlbg = w; w += ng_;
    m.wubg = w; w += ng_;
  }

  void KnitroInterface::solve(Memory& mem) const {
    KnitroMemory& m = dynamic_cast<KnitroMemory&>(mem);

    // Allocate KNITRO memory block (move back to init!)
    casadi_assert(m.kc_handle==0);
    m.kc_handle = KTR_new();
    casadi_assert(m.kc_handle!=0);
    int status;

    // Jacobian sparsity
    vector<int> Jcol, Jrow;
    if (!jacg_sp_.isNull()) {
      Jcol = jacg_sp_.get_col();
      Jrow = jacg_sp_.get_row();
    }

    // Hessian sparsity
    int nnzH = hesslag_sp_.isNull() ? 0 : hesslag_sp_.nnz();
    vector<int> Hcol, Hrow;
    if (nnzH>0) {
      Hcol = hesslag_sp_.get_col();
      Hrow = hesslag_sp_.get_row();
      status = KTR_set_int_param_by_name(m.kc_handle, "hessopt", KTR_HESSOPT_EXACT);
      casadi_assert_message(status==0, "KTR_set_int_param failed");
    } else {
      status = KTR_set_int_param_by_name(m.kc_handle, "hessopt", KTR_HESSOPT_LBFGS);
      casadi_assert_message(status==0, "KTR_set_int_param failed");
    }

    // Set user set options
    for (auto&& pp : double_param_) {
      status = KTR_set_double_param_by_name(m.kc_handle, pp.first.c_str(), pp.second);
      casadi_assert_message(status==0, "KnitroInterface::evaluate: cannot set " + pp.first);
    }

    for (auto&& pp : int_param_) {
      status = KTR_set_int_param_by_name(m.kc_handle, pp.first.c_str(), pp.second);
      casadi_assert_message(status==0, "KnitroInterface::evaluate: cannot set " + pp.first);
    }

    for (auto&& pp : string_param_) {
      status = KTR_set_char_param_by_name(m.kc_handle, pp.first.c_str(), pp.second.c_str());
      casadi_assert_message(status==0, "KnitroInterface::evaluate: cannot set " + pp.first);
    }

    // Type of constraints
    vector<int> cType(ng_, KTR_CONTYPE_GENERAL);
    if (hasSetOption("contype")) {
      vector<int> contype = option("contype");
      casadi_assert(contype.size()==cType.size());
      copy(contype.begin(), contype.end(), cType.begin());
    }

    // "Correct" upper and lower bounds
    casadi_copy(m.x0, nx_, m.wx);
    casadi_copy(m.lbx, nx_, m.wlbx);
    casadi_copy(m.ubx, nx_, m.wubx);
    casadi_copy(m.lbg, ng_, m.wlbg);
    casadi_copy(m.ubg, ng_, m.wubg);
    for (int i=0; i<nx_; ++i) if (isinf(m.wlbx[i])) m.wlbx[i] = -KTR_INFBOUND;
    for (int i=0; i<nx_; ++i) if (isinf(m.wubx[i])) m.wubx[i] =  KTR_INFBOUND;
    for (int i=0; i<ng_; ++i) if (isinf(m.wlbg[i])) m.wlbg[i] = -KTR_INFBOUND;
    for (int i=0; i<ng_; ++i) if (isinf(m.wubg[i])) m.wubg[i] =  KTR_INFBOUND;

    // Initialize KNITRO
    status = KTR_init_problem(m.kc_handle, nx_, KTR_OBJGOAL_MINIMIZE,
                              KTR_OBJTYPE_GENERAL, m.wlbx, m.wubx, ng_, get_ptr(cType),
                              m.wlbg, m.wubg, Jcol.size(), get_ptr(Jcol), get_ptr(Jrow),
                              nnzH, get_ptr(Hrow), get_ptr(Hcol), m.wx, 0); // initial lambda
    casadi_assert_message(status==0, "KTR_init_problem failed");

    // Register callback functions
    status = KTR_set_func_callback(m.kc_handle, &callback);
    casadi_assert_message(status==0, "KTR_set_func_callback failed");

    status = KTR_set_grad_callback(m.kc_handle, &callback);
    casadi_assert_message(status==0, "KTR_set_grad_callbackfailed");

    if (nnzH>0) {
      status = KTR_set_hess_callback(m.kc_handle, &callback);
      casadi_assert_message(status==0, "KTR_set_hess_callbackfailed");
    }

    // Lagrange multipliers
    vector<double> lambda(nx_+ng_);

    // Solve NLP
    double f;
    status = KTR_solve(m.kc_handle, m.wx, get_ptr(lambda), 0, &f,
                       0, 0, 0, 0, 0, static_cast<void*>(&m));
    m.return_status = return_codes(status);

    // Output primal solution
    casadi_copy(m.wx, nx_, m.x);

    // Output dual solution
    casadi_copy(get_ptr(lambda), ng_, m.lam_g);
    casadi_copy(get_ptr(lambda)+ng_, nx_, m.lam_x);

    // Output optimal cost
    if (m.f) *m.f = f;

    // Calculate constraints
    if (m.g) calc_fg(m, m.wx, m.p, 0, m.g);

    // Free memory (move to destructor!)
    KTR_free(&m.kc_handle);
    m.kc_handle = 0;
  }

  int KnitroInterface::callback(const int evalRequestCode, const int n, const int m, const int nnzJ,
                               const int nnzH, const double* const x, const double* const lambda,
                               double* const obj, double* const c, double* const objGrad,
                               double* const jac, double* const hessian, double* const hessVector,
                               void *userParams) {
    try {
      // Get a pointer to the calling object
      KnitroMemory& m = *static_cast<KnitroMemory*>(userParams);

      // Direct to the correct function
      switch (evalRequestCode) {
      case KTR_RC_EVALFC:
        m.self.calc_fg(m, x, m.p, obj, c);
        break;
      case KTR_RC_EVALGA:
        m.self.calc_gf_jg(m, x, m.p, objGrad, jac);
        break;
      case KTR_RC_EVALH:
        {
          double sigma = 1.;
          if (m.self.calc_hess_l(m, x, m.p, &sigma, lambda, hessian)) {
            casadi_error("calc_hess_l failed");
          }
        }
        break;
      default:
        casadi_error("KnitroInterface::callback: unknown method");
      }

      return 0;
    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "KnitroInterface::callback caugth exception: "
                               << ex.what() << endl;
      return -1;
    }
  }

  const char* KnitroInterface::return_codes(int flag) {
    switch (flag) {
    case KTR_RC_OPTIMAL_OR_SATISFACTORY: return "KTR_RC_OPTIMAL_OR_SATISFACTORY";
    case KTR_RC_NEAR_OPT: return "KTR_RC_NEAR_OPT";
    case KTR_RC_FEAS_XTOL: return "KTR_RC_FEAS_XTOL";
    case KTR_RC_FEAS_NO_IMPROVE: return "KTR_RC_FEAS_NO_IMPROVE";
    case KTR_RC_FEAS_FTOL: return "KTR_RC_FEAS_FTOL";
    case KTR_RC_INFEASIBLE: return "KTR_RC_INFEASIBLE";
    case KTR_RC_INFEAS_XTOL: return "KTR_RC_INFEAS_XTOL";
    case KTR_RC_INFEAS_NO_IMPROVE: return "KTR_RC_INFEAS_NO_IMPROVE";
    case KTR_RC_INFEAS_MULTISTART: return "KTR_RC_INFEAS_MULTISTART";
    case KTR_RC_INFEAS_CON_BOUNDS: return "KTR_RC_INFEAS_CON_BOUNDS";
    case KTR_RC_INFEAS_VAR_BOUNDS: return "KTR_RC_INFEAS_VAR_BOUNDS";
    case KTR_RC_UNBOUNDED: return "KTR_RC_UNBOUNDED";
    case KTR_RC_ITER_LIMIT_FEAS: return "KTR_RC_ITER_LIMIT_FEAS";
    case KTR_RC_TIME_LIMIT_FEAS: return "KTR_RC_TIME_LIMIT_FEAS";
    case KTR_RC_FEVAL_LIMIT_FEAS: return "KTR_RC_FEVAL_LIMIT_FEAS";
    case KTR_RC_MIP_EXH_FEAS: return "KTR_RC_MIP_EXH_FEAS";
    case KTR_RC_MIP_TERM_FEAS: return "KTR_RC_MIP_TERM_FEAS";
    case KTR_RC_MIP_SOLVE_LIMIT_FEAS: return "KTR_RC_MIP_SOLVE_LIMIT_FEAS";
    case KTR_RC_MIP_NODE_LIMIT_FEAS: return "KTR_RC_MIP_NODE_LIMIT_FEAS";
    case KTR_RC_ITER_LIMIT_INFEAS: return "KTR_RC_ITER_LIMIT_INFEAS";
    case KTR_RC_TIME_LIMIT_INFEAS: return "KTR_RC_TIME_LIMIT_INFEAS";
    case KTR_RC_FEVAL_LIMIT_INFEAS: return "KTR_RC_FEVAL_LIMIT_INFEAS";
    case KTR_RC_MIP_EXH_INFEAS: return "KTR_RC_MIP_EXH_INFEAS";
    case KTR_RC_MIP_SOLVE_LIMIT_INFEAS: return "KTR_RC_MIP_SOLVE_LIMIT_INFEAS";
    case KTR_RC_MIP_NODE_LIMIT_INFEAS: return "KTR_RC_MIP_NODE_LIMIT_INFEAS";
    case KTR_RC_CALLBACK_ERR: return "KTR_RC_CALLBACK_ERR";
    case KTR_RC_LP_SOLVER_ERR: return "KTR_RC_LP_SOLVER_ERR";
    case KTR_RC_EVAL_ERR: return "KTR_RC_EVAL_ERR";
    case KTR_RC_OUT_OF_MEMORY: return "KTR_RC_OUT_OF_MEMORY";
    case KTR_RC_USER_TERMINATION: return "KTR_RC_USER_TERMINATION";
    case KTR_RC_OPEN_FILE_ERR: return "KTR_RC_OPEN_FILE_ERR";
    case KTR_RC_BAD_N_OR_F: return "KTR_RC_BAD_N_OR_F";
    case KTR_RC_BAD_CONSTRAINT: return "KTR_RC_BAD_CONSTRAINT";
    case KTR_RC_BAD_JACOBIAN: return "KTR_RC_BAD_JACOBIAN";
    case KTR_RC_BAD_HESSIAN: return "KTR_RC_BAD_HESSIAN";
    case KTR_RC_BAD_CON_INDEX: return "KTR_RC_BAD_CON_INDEX";
    case KTR_RC_BAD_JAC_INDEX: return "KTR_RC_BAD_JAC_INDEX";
    case KTR_RC_BAD_HESS_INDEX: return "KTR_RC_BAD_HESS_INDEX";
    case KTR_RC_BAD_CON_BOUNDS: return "KTR_RC_BAD_CON_BOUNDS";
    case KTR_RC_BAD_VAR_BOUNDS: return "KTR_RC_BAD_VAR_BOUNDS";
    case KTR_RC_ILLEGAL_CALL: return "KTR_RC_ILLEGAL_CALL";
    case KTR_RC_BAD_KCPTR: return "KTR_RC_BAD_KCPTR";
    case KTR_RC_NULL_POINTER: return "KTR_RC_NULL_POINTER";
    case KTR_RC_BAD_INIT_VALUE: return "KTR_RC_BAD_INIT_VALUE";
    case KTR_RC_NEWPOINT_HALT: return "KTR_RC_NEWPOINT_HALT";
    case KTR_RC_BAD_LICENSE: return "KTR_RC_BAD_LICENSE";
    case KTR_RC_BAD_PARAMINPUT: return "KTR_RC_BAD_PARAMINPUT";
    case KTR_RC_LINEAR_SOLVER_ERR: return "KTR_RC_LINEAR_SOLVER_ERR";
    case KTR_RC_DERIV_CHECK_FAILED: return "KTR_RC_DERIV_CHECK_FAILED";
    case KTR_RC_DERIV_CHECK_TERMINATE: return "KTR_RC_DERIV_CHECK_TERMINATE";
    case KTR_RC_INTERNAL_ERROR: return "KTR_RC_INTERNAL_ERROR";
    }
    return 0;
  }

  KnitroMemory::KnitroMemory(const KnitroInterface& self) : self(self) {
    this->kc_handle = 0;
  }

  KnitroMemory::~KnitroMemory() {
    // Currently no persistent memory since KNITRO requires knowledge of nature of bounds
    // if (this->kc_handle) {
    //   KTR_free(&this->kc_handle);
    // }
  }

} // namespace casadi
