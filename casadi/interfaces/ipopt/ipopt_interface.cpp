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
#include "../../core/casadi_options.hpp"
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

    addOption("pass_nonlinear_variables", OT_BOOLEAN, false);
    addOption("print_time",               OT_BOOLEAN, true,
              "print information about execution time");

    // Monitors
    addOption("monitor",                  OT_STRINGVECTOR, GenericType(),  "",
              "eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h", true);

    // For passing metadata to IPOPT
    addOption("var_string_md",            OT_DICT, GenericType(),
              "String metadata (a dictionary with lists of strings) "
              "about variables to be passed to IPOPT");
    addOption("var_integer_md",           OT_DICT, GenericType(),
              "Integer metadata (a dictionary with lists of integers) "
              "about variables to be passed to IPOPT");
    addOption("var_numeric_md",           OT_DICT, GenericType(),
              "Numeric metadata (a dictionary with lists of reals) about "
              "variables to be passed to IPOPT");
    addOption("con_string_md",            OT_DICT, GenericType(),
              "String metadata (a dictionary with lists of strings) about "
              "constraints to be passed to IPOPT");
    addOption("con_integer_md",           OT_DICT, GenericType(),
              "Integer metadata (a dictionary with lists of integers) "
              "about constraints to be passed to IPOPT");
    addOption("con_numeric_md",           OT_DICT, GenericType(),
              "Numeric metadata (a dictionary with lists of reals) about "
              "constraints to be passed to IPOPT");

    // Start an application (temporarily)
    Ipopt::IpoptApplication temp_app;

    // Get all options available in (s)IPOPT
    map<string, Ipopt::SmartPtr<Ipopt::RegisteredOption> > regops =
        temp_app.RegOptions()->RegisteredOptionsList();
    for (map<string, Ipopt::SmartPtr<Ipopt::RegisteredOption> >::const_iterator it=regops.begin();
        it!=regops.end();
        ++it) {
      // Option identifier
      string opt_name = it->first;

      // Short description goes here, even though we do have a longer description
      string opt_desc = it->second->ShortDescription() + " (see IPOPT documentation)";

      // Get the type
      Ipopt::RegisteredOptionType ipopt_type = it->second->Type();
      TypeID casadi_type;

      // Map Ipopt option category to a CasADi options type
      switch (ipopt_type) {
      case Ipopt::OT_Number:    casadi_type = OT_REAL;          break;
      case Ipopt::OT_Integer:   casadi_type = OT_INTEGER;       break;
      case Ipopt::OT_String:    casadi_type = OT_STRING;        break;
      case Ipopt::OT_Unknown:   continue; // NOTE: No mechanism to handle OT_Unknown options
      default:                  continue; // NOTE: Unknown Ipopt options category
      }

      addOption(opt_name, casadi_type, GenericType(), opt_desc);

      // Set default values of IPOPT options
      if (casadi_type == OT_REAL) {
        setDefault(opt_name, it->second->DefaultNumber());
      } else if (casadi_type == OT_INTEGER) {
        setDefault(opt_name, it->second->DefaultInteger());
      } else if (casadi_type == OT_STRING) {
        setDefault(opt_name, it->second->DefaultString());
      }

      // Save to map containing IPOPT specific options
      ops_[opt_name] = casadi_type;
    }

    char * default_solver = getenv("IPOPT_DEFAULT_LINEAR_SOLVER");
    if (default_solver) {
      setOption("linear_solver", default_solver);
    }
  }

  IpoptInterface::~IpoptInterface() {
  }

  void IpoptInterface::init() {
    // Call the init method of the base class
    Nlpsol::init();

    // Read user options
    exact_hessian_ = !hasSetOption("hessian_approximation") ||
        option("hessian_approximation")=="exact";

    // Identify nonlinear variables
    pass_nonlinear_variables_ = option("pass_nonlinear_variables");

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

    bool ret = true;

    // Pass all the options to ipopt
    for (map<string, TypeID>::const_iterator it=ops_.begin(); it!=ops_.end(); ++it)
      if (hasSetOption(it->first)) {
        GenericType op = option(it->first);
        switch (it->second) {
        case OT_REAL:
          ret &= (*app)->Options()->SetNumericValue(it->first, op.toDouble(), false);
          break;
        case OT_INTEGER:
          ret &= (*app)->Options()->SetIntegerValue(it->first, op.toInt(), false);
          break;
        case OT_STRING:
          ret &= (*app)->Options()->SetStringValue(it->first, op.toString(), false);
          break;
        default:
          throw CasadiException("Illegal type");
        }
      }

    if (!ret) casadi_error("IpoptInterface::Init: Invalid options were detected by Ipopt.");

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
    m.t_mainloop = diffTimers(getTimerTime(), time0);

    // Save results to outputs
    casadi_copy(&m.fk, 1, m.f);
    casadi_copy(m.xk, nx_, m.x);
    casadi_copy(m.lam_gk, ng_, m.lam_g);
    casadi_copy(m.lam_xk, nx_, m.lam_x);
    casadi_copy(m.gk, ng_, m.g);

    if (status == Solve_Succeeded)
      m.return_status = "Solve_Succeeded";
    if (status == Solved_To_Acceptable_Level)
      m.return_status = "Solved_To_Acceptable_Level";
    if (status == Infeasible_Problem_Detected)
      m.return_status = "Infeasible_Problem_Detected";
    if (status == Search_Direction_Becomes_Too_Small)
      m.return_status = "Search_Direction_Becomes_Too_Small";
    if (status == Diverging_Iterates)
      m.return_status = "Diverging_Iterates";
    if (status == User_Requested_Stop)
      m.return_status = "User_Requested_Stop";
    if (status == Maximum_Iterations_Exceeded)
      m.return_status = "Maximum_Iterations_Exceeded";
    if (status == Restoration_Failed)
      m.return_status = "Restoration_Failed";
    if (status == Error_In_Step_Computation)
      m.return_status = "Error_In_Step_Computation";
    if (status == Not_Enough_Degrees_Of_Freedom)
      m.return_status = "Not_Enough_Degrees_Of_Freedom";
    if (status == Invalid_Problem_Definition)
      m.return_status = "Invalid_Problem_Definition";
    if (status == Invalid_Option)
      m.return_status = "Invalid_Option";
    if (status == Invalid_Number_Detected)
      m.return_status = "Invalid_Number_Detected";
    if (status == Unrecoverable_Exception)
      m.return_status = "Unrecoverable_Exception";
    if (status == NonIpopt_Exception_Thrown)
      m.return_status = "NonIpopt_Exception_Thrown";
    if (status == Insufficient_Memory)
      m.return_status = "Insufficient_Memory";
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
      if (!fcallback_.isNull()) {
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
              "You will only be able to use getStats(). "
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
      if (option("iteration_callback_ignore_errors")) {
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
    if (hasSetOption("var_string_md")) {
      Dict dict = option("var_string_md");
      for (Dict::const_iterator it=dict.begin(); it!=dict.end(); ++it) {
        string key = it->first; // Get the key
        vector<string> entry = it->second; // Get the entry
        // Check length for consistency
        casadi_assert_message(entry.size()==nx_, "Inconsistent length of IPOPT metadata.");
        var_string_md[key] = entry; // Save to IPOPT data structure
      }
    }

    if (hasSetOption("var_integer_md")) {
      Dict dict = option("var_integer_md");
      for (Dict::const_iterator it=dict.begin(); it!=dict.end(); ++it) {
        string key = it->first; // Get the key
        vector<int> entry = it->second; // Get the entry
        // Check length for consistency
        casadi_assert_message(entry.size()==nx_, "Inconsistent length of IPOPT metadata.");
        var_integer_md[key] = entry; // Save to IPOPT data structure
      }
    }

    if (hasSetOption("var_numeric_md")) {
      Dict dict = option("var_numeric_md");
      for (Dict::const_iterator it=dict.begin(); it!=dict.end(); ++it) {
        string key = it->first; // Get the key
        vector<double> entry = it->second; // Get the entry
        // Check length for consistency
        casadi_assert_message(entry.size()==nx_, "Inconsistent length of IPOPT metadata.");
        var_numeric_md[key] = entry; // Save to IPOPT data structure
      }
    }

    if (hasSetOption("con_string_md")) {
      Dict dict = option("con_string_md");
      for (Dict::const_iterator it=dict.begin(); it!=dict.end(); ++it) {
        string key = it->first; // Get the key
        vector<string> entry = it->second; // Get the entry
        // Check length for consistency
        casadi_assert_message(entry.size()==ng_, "Inconsistent length of IPOPT metadata.");
        con_string_md[key] = entry; // Save to IPOPT data structure
      }
    }

    if (hasSetOption("con_integer_md")) {
      Dict dict = option("con_integer_md");
      for (Dict::const_iterator it=dict.begin(); it!=dict.end(); ++it) {
        string key = it->first; // Get the key
        vector<int> entry = it->second; // Get the entry
        // Check length for consistency
        casadi_assert_message(entry.size()==ng_, "Inconsistent length of IPOPT metadata.");
        con_integer_md[key] = entry; // Save to IPOPT data structure
      }
    }

    if (hasSetOption("con_numeric_md")) {
      Dict dict = option("con_numeric_md");
      for (Dict::const_iterator it=dict.begin(); it!=dict.end(); ++it) {
        string key = it->first; // Get the key
        vector<double> entry = it->second; // Get the entry
        // Check length for consistency
        casadi_assert_message(entry.size()==ng_, "Inconsistent length of IPOPT metadata.");
        con_numeric_md[key] = entry; // Save to IPOPT data structure
      }
    }

    return true;
  }

  void IpoptInterface::setDefaultOptions(const std::vector<std::string>& recipes) {
    // Can be enabled when a new bugfixed version of Ipopt comes out
    //setOption("mehrotra_algorithm", "yes");
    //setOption("mu_oracle", "probing");
    for (int i=0;i<recipes.size();++i) {
      if (recipes[i]=="qp") {
        setOption("fixed_variable_treatment", "relax_bounds");
        setOption("jac_c_constant", "yes");
        setOption("jac_d_constant", "yes");
        setOption("hessian_constant", "yes");
      }
    }
  }

  IpoptMemory::IpoptMemory() {
    this->app = 0;
    this->userclass = 0;
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

} // namespace casadi
