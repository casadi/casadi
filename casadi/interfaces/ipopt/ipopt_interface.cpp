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



#include "ipopt_interface.hpp"
#include "ipopt_nlp.hpp"
#include "casadi/core/std_vector_tools.hpp"
#include "../../core/global_options.hpp"
#include "../../core/casadi_interrupt.hpp"

#include <ctime>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace std;
#include <IpIpoptApplication.hpp>

namespace casadi {
  extern "C"
  int CASADI_NLPSOL_IPOPT_EXPORT
  casadi_register_nlpsol_ipopt(Nlpsol::Plugin* plugin) {
    plugin->creator = IpoptInterface::creator;
    plugin->name = "ipopt";
    plugin->doc = IpoptInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_IPOPT_EXPORT casadi_load_nlpsol_ipopt() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_ipopt);
  }

  IpoptInterface::IpoptInterface(const std::string& name, const XProblem& nlp)
    : Nlpsol(name, nlp) {
  }

  IpoptInterface::~IpoptInterface() {
  }

  Options IpoptInterface::options_
  = {{&Nlpsol::options_},
     {{"pass_nonlinear_variables",
       {OT_BOOL,
        "Pass list of variables entering nonlinearly to IPOPT"}},
      {"print_time",
       {OT_BOOL,
        "print information about execution time"}},
      {"ipopt",
       {OT_DICT,
        "Options to be passed to IPOPT"}},
      {"var_string_md",
       {OT_DICT,
        "String metadata (a dictionary with lists of strings) "
        "about variables to be passed to IPOPT"}},
      {"var_integer_md",
       {OT_DICT,
        "Integer metadata (a dictionary with lists of integers) "
        "about variables to be passed to IPOPT"}},
      {"var_numeric_md",
       {OT_DICT,
        "Numeric metadata (a dictionary with lists of reals) about "
        "variables to be passed to IPOPT"}},
      {"con_string_md",
       {OT_DICT,
        "String metadata (a dictionary with lists of strings) about "
        "constraints to be passed to IPOPT"}},
      {"con_integer_md",
       {OT_DICT,
        "Integer metadata (a dictionary with lists of integers) "
        "about constraints to be passed to IPOPT"}},
      {"con_numeric_md",
       {OT_DICT,
        "Numeric metadata (a dictionary with lists of reals) about "
        "constraints to be passed to IPOPT"}}
     }
  };

  void IpoptInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    // Default options
    pass_nonlinear_variables_ = false;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="ipopt") {
        opts_ = op.second;
      } else if (op.first=="pass_nonlinear_variables") {
        pass_nonlinear_variables_ = op.second;
      } else if (op.first=="var_string_md") {
        var_string_md_ = op.second;
      } else if (op.first=="var_integer_md") {
        var_integer_md_ = op.second;
      } else if (op.first=="var_numeric_md") {
        var_numeric_md_ = op.second;
      } else if (op.first=="con_string_md") {
        con_string_md_ = op.second;
      } else if (op.first=="con_integer_md") {
        con_integer_md_ = op.second;
      } else if (op.first=="con_numeric_md") {
        con_numeric_md_ = op.second;
      }
    }

    // Do we need second order derivatives?
    exact_hessian_ = true;
    auto hessian_approximation = opts_.find("hessian_approximation");
    if (hessian_approximation!=opts_.end()) {
      exact_hessian_ = hessian_approximation->second == "exact";
    }

    // Setup NLP functions
    setup_f(); // Objective
    setup_g(); // Constraints
    setup_grad_f(); // Gradient of the objective
    setup_jac_g(); // Jacobian of the constraits
    if (exact_hessian_) {
      setup_hess_l(); // Hessian of the Lagrangian
    } else if (pass_nonlinear_variables_) {
      if (nlp2_.is_sx) {
        const Problem<SX>& nlp = nlp2_;
        SX fg = veccat(vector<SX>{nlp.out[NL_F], nlp.out[NL_G]});
        nl_ex_ = nl_var(fg, nlp.in[NL_X]);
      } else {
        const Problem<MX>& nlp = nlp2_;
        MX fg = veccat(vector<MX>{nlp.out[NL_F], nlp.out[NL_G]});
        nl_ex_ = nl_var(fg, nlp.in[NL_X]);
      }

    }

    // Allocate work vectors
    alloc_w(nx_, true); // xk_
    alloc_w(ng_, true); // lam_gk_
    alloc_w(nx_, true); // lam_xk_
    alloc_w(ng_, true); // gk_
    alloc_w(nx_, true); // grad_fk_
    alloc_w(jacg_sp_.nnz(), true); // jac_gk_
    if (exact_hessian_) {
      alloc_w(hesslag_sp_.nnz(), true); // hess_lk_
    }
  }

  void IpoptInterface::init_memory(Memory& mem) const {
    Nlpsol::init_memory(mem);
    IpoptMemory& m = dynamic_cast<IpoptMemory&>(mem);

    // Start an IPOPT application
    Ipopt::SmartPtr<Ipopt::IpoptApplication> *app = new Ipopt::SmartPtr<Ipopt::IpoptApplication>();
    m.app = static_cast<void*>(app);
    *app = new Ipopt::IpoptApplication(false);

    // Direct output through casadi::userOut()
    StreamJournal* jrnl_raw = new StreamJournal("console", J_ITERSUMMARY);
    jrnl_raw->SetOutputStream(&casadi::userOut());
    jrnl_raw->SetPrintLevel(J_DBG, J_NONE);
    SmartPtr<Journal> jrnl = jrnl_raw;
    (*app)->Jnlst()->AddJournal(jrnl);

    // Create an Ipopt user class -- need to use Ipopts spart pointer class
    Ipopt::SmartPtr<Ipopt::TNLP> *userclass = new Ipopt::SmartPtr<Ipopt::TNLP>();
    m.userclass = static_cast<void*>(userclass);
    *userclass = new IpoptUserClass(*this, m);

    if (verbose_) {
      userOut() << "There are " << nx_ << " variables and " << ng_ << " constraints." << endl;
      if (exact_hessian_) userOut() << "Using exact Hessian" << endl;
      else             userOut() << "Using limited memory Hessian approximation" << endl;
    }

    // Get all options available in (s)IPOPT
    auto regops = (*app)->RegOptions()->RegisteredOptionsList();

    // Pass all the options to ipopt
    for (auto&& op : opts_) {
      // Find the option
      auto regops_it = regops.find(op.first);
      if (regops_it==regops.end()) {
        casadi_error("No such IPOPT option: " + op.first);
      }

      // Get the type
      Ipopt::RegisteredOptionType ipopt_type = regops_it->second->Type();

      // Pass to IPOPT
      bool ret;
      switch (ipopt_type) {
      case Ipopt::OT_Number:
        ret = (*app)->Options()->SetNumericValue(op.first, op.second.to_double(), false);
        break;
      case Ipopt::OT_Integer:
        ret = (*app)->Options()->SetIntegerValue(op.first, op.second.to_int(), false);
        break;
      case Ipopt::OT_String:
        ret = (*app)->Options()->SetStringValue(op.first, op.second.to_string(), false);
        break;
      case Ipopt::OT_Unknown:
      default:
        casadi_warning("Cannot handle option \"" + op.first + "\", ignored");
        continue;
      }
      if (!ret) casadi_error("Invalid options were detected by Ipopt.");
    }

    // Override IPOPT's default linear solver
    if (opts_.find("linear_solver") == opts_.end()) {
      char * default_solver = getenv("IPOPT_DEFAULT_LINEAR_SOLVER");
      if (default_solver) {
        bool ret = (*app)->Options()->SetStringValue("linear_solver", default_solver, false);
        casadi_assert_message(ret, "Corrupted IPOPT_DEFAULT_LINEAR_SOLVER environmental variable");
      }
    }

    // Intialize the IpoptApplication and process the options
    Ipopt::ApplicationReturnStatus status = (*app)->Initialize();
    casadi_assert_message(status == Solve_Succeeded, "Error during IPOPT initialization");
  }

  void IpoptInterface::set_work(Memory& mem, const double**& arg, double**& res,
                                int*& iw, double*& w) const {
    IpoptMemory& m = dynamic_cast<IpoptMemory&>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    // Work vectors
    m.xk = w; w += nx_;
    m.lam_gk = w; w += ng_;
    m.lam_xk = w; w += nx_;
    m.gk = w; w += ng_;
    m.grad_fk = w; w += nx_;
    m.jac_gk = w; w += jacg_sp_.nnz();
    if (exact_hessian_) {
      m.hess_lk = w; w += hesslag_sp_.nnz();
    }
  }

  inline const char* return_status_string(Ipopt::ApplicationReturnStatus status) {
    switch (status) {
    case Solve_Succeeded:
      return "Solve_Succeeded";
    case Solved_To_Acceptable_Level:
      return "Solved_To_Acceptable_Level";
    case Infeasible_Problem_Detected:
      return "Infeasible_Problem_Detected";
    case Search_Direction_Becomes_Too_Small:
      return "Search_Direction_Becomes_Too_Small";
    case Diverging_Iterates:
      return "Diverging_Iterates";
    case User_Requested_Stop:
      return "User_Requested_Stop";
    case Maximum_Iterations_Exceeded:
      return "Maximum_Iterations_Exceeded";
    case Restoration_Failed:
      return "Restoration_Failed";
    case Error_In_Step_Computation:
      return "Error_In_Step_Computation";
    case Not_Enough_Degrees_Of_Freedom:
      return "Not_Enough_Degrees_Of_Freedom";
    case Invalid_Problem_Definition:
      return "Invalid_Problem_Definition";
    case Invalid_Option:
      return "Invalid_Option";
    case Invalid_Number_Detected:
      return "Invalid_Number_Detected";
    case Unrecoverable_Exception:
      return "Unrecoverable_Exception";
    case NonIpopt_Exception_Thrown:
      return "NonIpopt_Exception_Thrown";
    case Insufficient_Memory:
      return "Insufficient_Memory";
    case Internal_Error:
      return "Internal_Error";
    case Maximum_CpuTime_Exceeded:
      return "Maximum_CpuTime_Exceeded";
    case Feasible_Point_Found:
      return "Feasible_Point_Found";
    }
    return "Unknown";
  }

  void IpoptInterface::solve(Memory& mem) const {
    IpoptMemory& m = dynamic_cast<IpoptMemory&>(mem);

    // Check the provided inputs
    checkInputs(mem);

    // Reset statistics
    if (gather_stats_) {
      m.inf_pr.clear();
      m.inf_du.clear();
      m.mu.clear();
      m.d_norm.clear();
      m.regularization_size.clear();
      m.alpha_pr.clear();
      m.alpha_du.clear();
      m.obj.clear();
      m.ls_trials.clear();
    }

    // Reset number of iterations
    m.n_iter = 0;

    // Reset function timers
    m.t_calc_f = m.t_calc_g = m.t_calc_grad_f = m.t_calc_jac_g = m.t_calc_hess_l = 0;

    // Reset function counters
    m.n_calc_f = m.n_calc_g = m.n_calc_grad_f = m.n_calc_jac_g = m.n_calc_hess_l = 0;

    // Legacy
    m.t_callback_fun = m.t_callback_prepare = m.t_mainloop = {0, 0};
    m.n_eval_callback = 0;

    // Get back the smart pointers
    Ipopt::SmartPtr<Ipopt::TNLP> *userclass =
      static_cast<Ipopt::SmartPtr<Ipopt::TNLP>*>(m.userclass);
    Ipopt::SmartPtr<Ipopt::IpoptApplication> *app =
      static_cast<Ipopt::SmartPtr<Ipopt::IpoptApplication>*>(m.app);

    Timer time0 = getTimerTime();
    // Ask Ipopt to solve the problem
    Ipopt::ApplicationReturnStatus status = (*app)->OptimizeTNLP(*userclass);
    m.return_status = return_status_string(status);
    m.t_mainloop = diffTimers(getTimerTime(), time0);

    // Save results to outputs
    casadi_copy(&m.fk, 1, m.f);
    casadi_copy(m.xk, nx_, m.x);
    casadi_copy(m.lam_gk, ng_, m.lam_g);
    casadi_copy(m.lam_xk, nx_, m.lam_x);
    casadi_copy(m.gk, ng_, m.g);
  }

  bool IpoptInterface::
  intermediate_callback(IpoptMemory& m, const double* x, const double* z_L, const double* z_U,
                        const double* g, const double* lambda, double obj_value, int iter,
                        double inf_pr, double inf_du, double mu, double d_norm,
                        double regularization_size, double alpha_du, double alpha_pr,
                        int ls_trials, bool full_callback) const {
    m.n_iter += 1;
    try {
      log("intermediate_callback started");
      if (gather_stats_) {
        m.inf_pr.push_back(inf_pr);
        m.inf_du.push_back(inf_du);
        m.mu.push_back(mu);
        m.d_norm.push_back(d_norm);
        m.regularization_size.push_back(regularization_size);
        m.alpha_pr.push_back(alpha_pr);
        m.alpha_du.push_back(alpha_du);
        m.ls_trials.push_back(ls_trials);
        m.obj.push_back(obj_value);
      }
      Timer time0 = getTimerTime();
      if (!fcallback_.is_null()) {
        if (full_callback) {
          casadi_copy(x, nx_, m.xk);
          for (int i=0; i<nx_; ++i) {
            m.lam_xk[i] = z_U[i]-z_L[i];
          }
          casadi_copy(lambda, ng_, m.lam_gk);
          casadi_copy(g, ng_, m.gk);
        } else {
          if (iter==0) {
            userOut<true, PL_WARN>()
              << "Warning: intermediate_callback is disfunctional in your installation. "
              "You will only be able to use stats(). "
              "See https://github.com/casadi/casadi/wiki/enableIpoptCallback to enable it."
              << endl;
          }
        }

        // Inputs
        fill_n(m.arg, fcallback_.n_in(), nullptr);
        m.arg[NLPSOL_X] = x;
        m.arg[NLPSOL_F] = &obj_value;
        m.arg[NLPSOL_G] = g;
        m.arg[NLPSOL_LAM_P] = 0;
        m.arg[NLPSOL_LAM_X] = m.lam_xk;
        m.arg[NLPSOL_LAM_G] = m.lam_gk;

        // Outputs
        fill_n(m.res, fcallback_.n_out(), nullptr);
        double ret_double;
        m.res[0] = &ret_double;

        // Evaluate the callback function
        m.n_eval_callback += 1;
        fcallback_(m.arg, m.res, m.iw, m.w, 0);
        int ret = static_cast<int>(ret_double);

        DiffTime delta = diffTimers(getTimerTime(), time0);
        timerPlusEq(m.t_callback_fun, delta);
        return  !ret;
      } else {
        return 1;
      }
    } catch(exception& ex) {
      if (iteration_callback_ignore_errors_) {
        userOut<true, PL_WARN>() << "intermediate_callback: " << ex.what() << endl;
      } else {
        throw ex;
      }
      return 1;
    }
  }

  void IpoptInterface::
  finalize_solution(IpoptMemory& m, const double* x, const double* z_L, const double* z_U,
                    const double* g, const double* lambda, double obj_value,
                    int iter_count) const {
    try {
      // Get primal solution
      casadi_copy(x, nx_, m.xk);

      // Get optimal cost
      m.fk = obj_value;

      // Get dual solution (simple bounds)
      if (m.lam_xk) {
        for (int i=0; i<nx_; ++i) {
          m.lam_xk[i] = z_U[i]-z_L[i];
        }
      }

      // Get dual solution (nonlinear bounds)
      casadi_copy(lambda, ng_, m.lam_gk);

      // Get the constraints
      casadi_copy(g, ng_, m.gk);

      // Get statistics
      m.iter_count = iter_count;

    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "finalize_solution failed: " << ex.what() << endl;
    }
  }

  bool IpoptInterface::
  get_bounds_info(IpoptMemory& m, double* x_l, double* x_u,
                  double* g_l, double* g_u) const {
    try {
      casadi_copy(m.lbx, nx_, x_l);
      casadi_copy(m.ubx, nx_, x_u);
      casadi_copy(m.lbg, ng_, g_l);
      casadi_copy(m.ubg, ng_, g_u);
      return true;
    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "get_bounds_info failed: " << ex.what() << endl;
      return false;
    }
  }

  bool IpoptInterface::
  get_starting_point(IpoptMemory& m, bool init_x, double* x,
                     bool init_z, double* z_L, double* z_U,
                     bool init_lambda, double* lambda) const {
    try {
      // Initialize primal variables
      if (init_x) {
        casadi_copy(m.x0, nx_, x);
      }

      // Initialize dual variables (simple bounds)
      if (init_z) {
        if (m.lam_x0) {
          for (int i=0; i<nx_; ++i) {
            z_L[i] = max(0., -m.lam_x0[i]);
            z_U[i] = max(0., m.lam_x0[i]);
          }
        } else {
          casadi_fill(z_L, nx_, 0.);
          casadi_fill(z_U, nx_, 0.);
        }
      }

      // Initialize dual variables (nonlinear bounds)
      if (init_lambda) {
        casadi_copy(m.lam_g0, ng_, lambda);
      }

      return true;
    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "get_starting_point failed: " << ex.what() << endl;
      return false;
    }
  }

  void IpoptInterface::get_nlp_info(IpoptMemory& m, int& nx, int& ng,
                                    int& nnz_jac_g, int& nnz_h_lag) const {
    try {
      // Number of variables
      nx = nx_;

      // Number of constraints
      ng = ng_;

      // Number of Jacobian nonzeros
      nnz_jac_g = ng_==0 ? 0 : jacg_sp_.nnz();

      // Number of Hessian nonzeros (only upper triangular half)
      nnz_h_lag = exact_hessian_ ? hesslag_sp_.nnz() : 0;

    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "get_nlp_info failed: " << ex.what() << endl;
    }
  }

  int IpoptInterface::get_number_of_nonlinear_variables() const {
    try {
      if (!pass_nonlinear_variables_) {
        // No Hessian has been interfaced
        return -1;
      } else {
        // Number of variables that appear nonlinearily
        int nv = 0;
        for (auto&& i : nl_ex_) if (i) nv++;
        return nv;
      }
    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "get_number_of_nonlinear_variables failed: " << ex.what() << endl;
      return -1;
    }
  }

  bool IpoptInterface::
  get_list_of_nonlinear_variables(int num_nonlin_vars, int* pos_nonlin_vars) const {
    try {
      for (int i=0; i<nl_ex_.size(); ++i) {
        if (nl_ex_[i]) *pos_nonlin_vars++ = i;
      }
      return true;
    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "get_list_of_nonlinear_variables failed: " << ex.what() << endl;
      return false;
    }
  }

  bool IpoptInterface::
  get_var_con_metadata(map<string, vector<string> >& var_string_md,
                       map<string, vector<int> >& var_integer_md,
                       map<string, vector<double> >& var_numeric_md,
                       map<string, vector<string> >& con_string_md,
                       map<string, vector<int> >& con_integer_md,
                       map<string, vector<double> >& con_numeric_md) const {
    for (auto&& op : var_string_md_) var_string_md[op.first] = op.second;
    for (auto&& op : var_integer_md_) var_integer_md[op.first] = op.second;
    for (auto&& op : var_numeric_md_) var_numeric_md[op.first] = op.second;
    for (auto&& op : con_string_md_) con_string_md[op.first] = op.second;
    for (auto&& op : con_integer_md_) con_integer_md[op.first] = op.second;
    for (auto&& op : con_numeric_md_) con_numeric_md[op.first] = op.second;
    return true;
  }

  IpoptMemory::IpoptMemory() {
    this->app = 0;
    this->userclass = 0;
    this->return_status = "Unset";
  }

  IpoptMemory::~IpoptMemory() {
    // Free Ipopt application instance (or rather, the smart pointer holding it)
    if (this->app != 0) {
      delete static_cast<Ipopt::SmartPtr<Ipopt::IpoptApplication>*>(this->app);
    }

    // Free Ipopt user class (or rather, the smart pointer holding it)
    if (this->userclass != 0) {
      delete static_cast<Ipopt::SmartPtr<Ipopt::TNLP>*>(this->userclass);
    }
  }

  Dict IpoptMemory::get_stats() const {
    Dict stats = NlpsolMemory::get_stats();
    stats["return_status"] = return_status;
    stats["iter_count"] = iter_count;
    return stats;
  }

} // namespace casadi
