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

    addOption("knitro", OT_DICT,
              "Options to be passed to KNITRO");
    addOption("contype", OT_INTVECTOR,
              "Type of constraint");
  }


  KnitroInterface::~KnitroInterface() {
  }

  void KnitroInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="knitro") {
        opts_ = op.second;
      } else if (op.first=="contype") {
        contype_ = op.second;
      }
    }

    // Type of constraints, general by default
    if (contype_.empty()) {
      contype_.resize(ng_, KTR_CONTYPE_GENERAL);
    } else {
      casadi_assert(contype_.size()==ng_);
    }

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

    // Pass user set options
    for (auto&& op : opts_) {
      // Try double
      if (op.second.can_cast_to(OT_DOUBLE)) {
        status = KTR_set_double_param_by_name(m.kc_handle, op.first.c_str(), op.second);
        if (status==0) continue;
      }

      // Try integer
      if (op.second.can_cast_to(OT_INT)) {
        status = KTR_set_int_param_by_name(m.kc_handle, op.first.c_str(), op.second);
        if (status==0) continue;
      }

      // try string
      if (op.second.can_cast_to(OT_STRING)) {
        string str = op.second.to_string();
        status = KTR_set_char_param_by_name(m.kc_handle, op.first.c_str(), str.c_str());
        if (status==0) continue;
      }

      // Error if reached this point
      casadi_error("KNITRO error setting option \"" + op.first + "\"");
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
                              KTR_OBJTYPE_GENERAL, m.wlbx, m.wubx, ng_, get_ptr(contype_),
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
