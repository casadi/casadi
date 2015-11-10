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

// Headers for sIPOPT
#ifdef WITH_SIPOPT
#include <SensApplication.hpp>
#include <IpPDSearchDirCalc.hpp>
#include <IpIpoptAlg.hpp>
#include <SensRegOp.hpp>
#include <SensReducedHessianCalculator.hpp>
#endif // WITH_SIPOPT

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

    // Set pointers to zero
    app_ = 0;
    userclass_ = 0;
#ifdef WITH_SIPOPT
    app_sens_ = 0;
#endif // WITH_SIPOPT

    // Start an application (temporarily)
    Ipopt::IpoptApplication temp_app;

    // Start a sensitivity application (temporarily)
#ifdef WITH_SIPOPT
    Ipopt::SensApplication temp_sens_app(temp_app.Jnlst(), temp_app.Options(),
                                         temp_app.RegOptions());

    // Register sIPOPT options
    Ipopt::RegisterOptions_sIPOPT(temp_app.RegOptions());
    temp_app.Options()->SetRegisteredOptions(temp_app.RegOptions());
#endif // WITH_SIPOPT

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

  void IpoptInterface::freeIpopt() {
    // Free sensitivity application (or rather, the smart pointer holding it)
#ifdef WITH_SIPOPT
    if (app_sens_ != 0) {
      delete static_cast<Ipopt::SmartPtr<Ipopt::SensApplication>*>(app_sens_);
      app_sens_ = 0;
    }
#endif // WITH_SIPOPT

    // Free Ipopt application instance (or rather, the smart pointer holding it)
    if (app_ != 0) {
      delete static_cast<Ipopt::SmartPtr<Ipopt::IpoptApplication>*>(app_);
      app_ = 0;
    }

    // Free Ipopt user class (or rather, the smart pointer holding it)
    if (userclass_ != 0) {
      delete static_cast<Ipopt::SmartPtr<Ipopt::TNLP>*>(userclass_);
      userclass_ = 0;
    }
  }

  IpoptInterface::~IpoptInterface() {
    freeIpopt();
  }

  void IpoptInterface::init() {
    // Free existing IPOPT instance
    freeIpopt();

    // Call the init method of the base class
    Nlpsol::init();

    // Read user options
    exact_hessian_ = !hasSetOption("hessian_approximation") ||
        option("hessian_approximation")=="exact";
#ifdef WITH_SIPOPT
    if (hasSetOption("run_sens")) {
      run_sens_ = option("run_sens")=="yes";
    } else {
      run_sens_  = false;
    }
    if (hasSetOption("compute_red_hessian")) {
      compute_red_hessian_ = option("compute_red_hessian")=="yes";
    } else {
      compute_red_hessian_ = false;
    }
#endif // WITH_SIPOPT

    // Get/generate required functions
    gradF();
    jacG();
    if (exact_hessian_) {
      hessLag();
    }

    // Start an IPOPT application
    Ipopt::SmartPtr<Ipopt::IpoptApplication> *app = new Ipopt::SmartPtr<Ipopt::IpoptApplication>();
    app_ = static_cast<void*>(app);
    *app = new Ipopt::IpoptApplication(false);

    // Direct output through casadi::userOut()
    StreamJournal* jrnl_raw = new StreamJournal("console", J_ITERSUMMARY);
    jrnl_raw->SetOutputStream(&casadi::userOut());
    jrnl_raw->SetPrintLevel(J_DBG, J_NONE);
    SmartPtr<Journal> jrnl = jrnl_raw;
    (*app)->Jnlst()->AddJournal(jrnl);

#ifdef WITH_SIPOPT
    if (run_sens_ || compute_red_hessian_) {
      // Start an sIPOPT application
      Ipopt::SmartPtr<Ipopt::SensApplication> *app_sens =
          new Ipopt::SmartPtr<Ipopt::SensApplication>();
      app_sens_ = static_cast<void*>(app_sens);
      *app_sens =
          new Ipopt::SensApplication((*app)->Jnlst(), (*app)->Options(), (*app)->RegOptions());

      // Register sIPOPT options
      Ipopt::RegisterOptions_sIPOPT((*app)->RegOptions());
      (*app)->Options()->SetRegisteredOptions((*app)->RegOptions());
    }
#endif // WITH_SIPOPT

    // Create an Ipopt user class -- need to use Ipopts spart pointer class
    Ipopt::SmartPtr<Ipopt::TNLP> *userclass = new Ipopt::SmartPtr<Ipopt::TNLP>();
    userclass_ = static_cast<void*>(userclass);
    *userclass = new IpoptUserClass(this);

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

    // Extra initialization required by sIPOPT
    //   #ifdef WITH_SIPOPT
    //   if (run_sens_ || compute_red_hessian_) {
    //     Ipopt::ApplicationReturnStatus status = (*app)->Initialize("");
    //     casadi_assert_message(status == Solve_Succeeded, "Error during IPOPT initialization");
    //   }
    //   #endif // WITH_SIPOPT

    // Intialize the IpoptApplication and process the options
    Ipopt::ApplicationReturnStatus status = (*app)->Initialize();
    casadi_assert_message(status == Solve_Succeeded, "Error during IPOPT initialization");

#ifdef WITH_SIPOPT
    if (run_sens_ || compute_red_hessian_) {
      Ipopt::SmartPtr<Ipopt::SensApplication> *app_sens =
          static_cast<Ipopt::SmartPtr<Ipopt::SensApplication> *>(app_sens_);
      (*app_sens)->Initialize();
    }
#endif // WITH_SIPOPT

    // Setup NLP functions
    if (nlp2_.is_sx) {
      setup<SX>();
    } else {
      setup<MX>();
    }

    // Allocate work vectors
    alloc_w(nx_, true); // xk_
    alloc_w(ng_, true); // lam_gk_
    alloc_w(ng_, true); // gk_
    alloc_w(nx_, true); // grad_fk_
    alloc_w(jacg_sp_.nnz(), true); // jac_gk_
    if (exact_hessian_) {
      alloc_w(hesslag_sp_.nnz(), true); // hess_lk_
    }
  }

  void IpoptInterface::reset(void* mem, const double**& arg, double**& res, int*& iw, double*& w) {
    // Reset the base classes
    Nlpsol::reset(mem, arg, res, iw, w);

    // Work vectors
    xk_ = w; w += nx_;
    lam_gk_ = w; w += ng_;
    gk_ = w; w += ng_;
    grad_fk_ = w; w += nx_;
    jac_gk_ = w; w += jacg_sp_.nnz();
    if (exact_hessian_) {
      hess_lk_ = w; w += hesslag_sp_.nnz();
    }

    // New iterate
    new_x_ = new_lam_f_ = new_lam_g_ = true;
  }

  void IpoptInterface::solve(void* mem) {
    for (int i=0; i<NLPSOL_NUM_IN; ++i) {
      const double *v;
      switch (i) {
      case NLPSOL_X0: v = x0_; break;
      case NLPSOL_P: v = p_; break;
      case NLPSOL_LBX: v = lbx_; break;
      case NLPSOL_UBX: v = ubx_; break;
      case NLPSOL_LBG: v = lbg_; break;
      case NLPSOL_UBG: v = ubg_; break;
      case NLPSOL_LAM_X0: v = lam_x0_; break;
      case NLPSOL_LAM_G0: v = lam_g0_; break;
      default: casadi_assert(0);
      }
      if (v) {
        setInputNZ(v, i);
      } else {
        setInput(0., i);
      }
    }

    // Check the provided inputs
    checkInputs(mem);

    if (gather_stats_) {
      Dict iterations;
      iterations["inf_pr"] = std::vector<double>();
      iterations["inf_du"] = std::vector<double>();
      iterations["mu"] = std::vector<double>();
      iterations["d_norm"] = std::vector<double>();
      iterations["regularization_size"] = std::vector<double>();
      iterations["obj"] = std::vector<double>();
      iterations["ls_trials"] = std::vector<int>();
      iterations["alpha_pr"] = std::vector<double>();
      iterations["alpha_du"] = std::vector<double>();
      iterations["obj"] = std::vector<double>();
      stats_["iterations"] = iterations;
    }

    // Reset number of iterations
    n_iter_ = 0;

    // Reset function timers
    t_calc_f_ = t_calc_g_ = t_calc_grad_f_ = t_calc_jac_g_ = t_calc_hess_l_ = 0;

    // Reset function counters
    n_calc_f_ = n_calc_g_ = n_calc_grad_f_ = n_calc_jac_g_ = n_calc_hess_l_ = 0;

    // Legacy
    t_callback_fun_ = t_callback_prepare_ = t_mainloop_ = {0, 0};
    n_eval_callback_ = 0;

    // Get back the smart pointers
    Ipopt::SmartPtr<Ipopt::TNLP> *userclass =
        static_cast<Ipopt::SmartPtr<Ipopt::TNLP>*>(userclass_);
    Ipopt::SmartPtr<Ipopt::IpoptApplication> *app =
        static_cast<Ipopt::SmartPtr<Ipopt::IpoptApplication>*>(app_);

    Timer time0 = getTimerTime();
    // Ask Ipopt to solve the problem
    Ipopt::ApplicationReturnStatus status = (*app)->OptimizeTNLP(*userclass);
    t_mainloop_ = diffTimers(getTimerTime(), time0);

#ifdef WITH_SIPOPT
    if (run_sens_ || compute_red_hessian_) {
      // Calculate parametric sensitivities
      Ipopt::SmartPtr<Ipopt::SensApplication> *app_sens =
          static_cast<Ipopt::SmartPtr<Ipopt::SensApplication>*>(app_sens_);
      (*app_sens)->SetIpoptAlgorithmObjects(*app, status);
      (*app_sens)->Run();

      // Access the reduced Hessian calculator
#ifdef WITH_CASADI_PATCH
      if (compute_red_hessian_) {
        // Get the reduced Hessian
        std::vector<double> red_hess = (*app_sens)->ReducedHessian();

        // Get the dimensions
        int N;
        for (N=0; N*N<red_hess.size(); ++N) {}
        casadi_assert(N*N==red_hess.size());

        // Store to statistics
        red_hess_ = DMatrix(Sparsity::dense(N, N), red_hess);
      }
#endif // WITH_CASADI_PATCH
    }
#endif // WITH_SIPOPT

    if (status == Solve_Succeeded)
      stats_["return_status"] = "Solve_Succeeded";
    if (status == Solved_To_Acceptable_Level)
      stats_["return_status"] = "Solved_To_Acceptable_Level";
    if (status == Infeasible_Problem_Detected)
      stats_["return_status"] = "Infeasible_Problem_Detected";
    if (status == Search_Direction_Becomes_Too_Small)
      stats_["return_status"] = "Search_Direction_Becomes_Too_Small";
    if (status == Diverging_Iterates)
      stats_["return_status"] = "Diverging_Iterates";
    if (status == User_Requested_Stop)
      stats_["return_status"] = "User_Requested_Stop";
    if (status == Maximum_Iterations_Exceeded)
      stats_["return_status"] = "Maximum_Iterations_Exceeded";
    if (status == Restoration_Failed)
      stats_["return_status"] = "Restoration_Failed";
    if (status == Error_In_Step_Computation)
      stats_["return_status"] = "Error_In_Step_Computation";
    if (status == Not_Enough_Degrees_Of_Freedom)
      stats_["return_status"] = "Not_Enough_Degrees_Of_Freedom";
    if (status == Invalid_Problem_Definition)
      stats_["return_status"] = "Invalid_Problem_Definition";
    if (status == Invalid_Option)
      stats_["return_status"] = "Invalid_Option";
    if (status == Invalid_Number_Detected)
      stats_["return_status"] = "Invalid_Number_Detected";
    if (status == Unrecoverable_Exception)
      stats_["return_status"] = "Unrecoverable_Exception";
    if (status == NonIpopt_Exception_Thrown)
      stats_["return_status"] = "NonIpopt_Exception_Thrown";
    if (status == Insufficient_Memory)
      stats_["return_status"] = "Insufficient_Memory";

    stats_["n_calc_f"] = n_calc_f_;
    stats_["t_calc_f"] = t_calc_f_;
    stats_["n_calc_g"] = n_calc_g_;
    stats_["t_calc_g"] = t_calc_g_;
    stats_["n_calc_grad_f"] = n_calc_grad_f_;
    stats_["t_calc_grad_f"] = t_calc_grad_f_;
    stats_["n_calc_jac_g"] = n_calc_jac_g_;
    stats_["t_calc_jac_g"] = t_calc_jac_g_;
    stats_["n_calc_hess_l"] = n_calc_hess_l_;
    stats_["t_calc_hess_l"] = t_calc_hess_l_;

    stats_["t_mainloop"] = diffToDict(t_mainloop_);
    stats_["t_callback_fun"] = diffToDict(t_callback_fun_);
    stats_["t_callback_prepare"] = diffToDict(t_callback_prepare_);
    stats_["n_eval_callback"] = n_eval_callback_;

    stats_["iter_count"] = n_iter_-1;

    for (int i=0; i<NLPSOL_NUM_OUT; ++i) {
      double **v;
      switch (i) {
      case NLPSOL_X: v = &x_; break;
      case NLPSOL_F: v = &f_; break;
      case NLPSOL_G: v = &g_; break;
      case NLPSOL_LAM_X: v = &lam_x_; break;
      case NLPSOL_LAM_G: v = &lam_g_; break;
      case NLPSOL_LAM_P: v = &lam_p_; break;
      default: casadi_assert(0);
      }
      if (*v) getOutputNZ(*v, i);
    }
  }

  bool IpoptInterface::intermediate_callback(
      const double* x, const double* z_L, const double* z_U, const double* g,
      const double* lambda, double obj_value, int iter,
      double inf_pr, double inf_du, double mu, double d_norm,
      double regularization_size, double alpha_du, double alpha_pr, int ls_trials,
      bool full_callback) {
    n_iter_ += 1;
    try {
      log("intermediate_callback started");
      if (gather_stats_) {
        Dict iterations = stats_["iterations"];
        append_to_vec(iterations["inf_pr"], inf_pr);
        append_to_vec(iterations["inf_du"], inf_du);
        append_to_vec(iterations["mu"], mu);
        append_to_vec(iterations["d_norm"], d_norm);
        append_to_vec(iterations["regularization_size"], regularization_size);
        append_to_vec(iterations["alpha_pr"], alpha_pr);
        append_to_vec(iterations["alpha_du"], alpha_du);
        append_to_vec(iterations["ls_trials"], ls_trials);
        append_to_vec(iterations["obj"], obj_value);
        stats_["iterations"] = iterations;
      }
      Timer time0 = getTimerTime();
      if (!fcallback_.isNull()) {
        if (full_callback) {
          if (!output(NLPSOL_X).is_empty()) copy(x, x+nx_, output(NLPSOL_X)->begin());

          vector<double>& lambda_x = output(NLPSOL_LAM_X).data();
          for (int i=0; i<lambda_x.size(); ++i) {
            lambda_x[i] = z_U[i]-z_L[i];
          }
          if (!output(NLPSOL_LAM_G).is_empty())
            copy(lambda, lambda+ng_, output(NLPSOL_LAM_G)->begin());
          if (!output(NLPSOL_G).is_empty()) copy(g, g+ng_, output(NLPSOL_G)->begin());
        } else {
          if (iter==0) {
            userOut<true, PL_WARN>()
              << "Warning: intermediate_callback is disfunctional in your installation. "
              "You will only be able to use getStats(). "
              "See https://github.com/casadi/casadi/wiki/enableIpoptCallback to enable it."
              << endl;
          }
        }

        Dict iteration;
        iteration["iter"] = iter;
        iteration["inf_pr"] = inf_pr;
        iteration["inf_du"] = inf_du;
        iteration["mu"] = mu;
        iteration["d_norm"] = d_norm;
        iteration["regularization_size"] = regularization_size;
        iteration["alpha_pr"] = alpha_pr;
        iteration["alpha_du"] = alpha_du;
        iteration["ls_trials"] = ls_trials;
        iteration["obj"] = obj_value;
        stats_["iteration"] = iteration;

        // Inputs
        fill_n(arg_, fcallback_.n_in(), nullptr);
        arg_[NLPSOL_X] = x;
        arg_[NLPSOL_F] = &obj_value;
        arg_[NLPSOL_G] = g;
        arg_[NLPSOL_LAM_P] = output(NLPSOL_LAM_P).ptr();
        arg_[NLPSOL_LAM_X] = output(NLPSOL_LAM_X).ptr();
        arg_[NLPSOL_LAM_G] = output(NLPSOL_LAM_G).ptr();

        // Outputs
        fill_n(res_, fcallback_.n_out(), nullptr);
        double ret_double;
        res_[0] = &ret_double;

        // Evaluate the callback function
        n_eval_callback_ += 1;
        fcallback_(arg_, res_, iw_, w_, 0);
        int ret = static_cast<int>(ret_double);

        DiffTime delta = diffTimers(getTimerTime(), time0);
        timerPlusEq(t_callback_fun_, delta);
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

  void IpoptInterface::finalize_solution(const double* x, const double* z_L, const double* z_U,
                                        const double* g, const double* lambda, double obj_value,
                                        int iter_count) {
    try {
      // Get primal solution
      copy(x, x+nx_, output(NLPSOL_X)->begin());

      // Get optimal cost
      output(NLPSOL_F).at(0) = obj_value;

      // Get dual solution (simple bounds)
      vector<double>& lambda_x = output(NLPSOL_LAM_X).data();
      for (int i=0; i<lambda_x.size(); ++i) {
        lambda_x[i] = z_U[i]-z_L[i];
      }

      // Get dual solution (nonlinear bounds)
      copy(lambda, lambda+ng_, output(NLPSOL_LAM_G)->begin());

      // Get the constraints
      copy(g, g+ng_, output(NLPSOL_G)->begin());

      // Get statistics
      stats_["iter_count"] = iter_count;

    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "finalize_solution failed: " << ex.what() << endl;
    }
  }

  bool IpoptInterface::get_bounds_info(int n, double* x_l, double* x_u,
                                      int m, double* g_l, double* g_u) {
    try {
      casadi_assert(n == nx_);
      casadi_assert(m == ng_);
      input(NLPSOL_LBX).getNZ(x_l);
      input(NLPSOL_UBX).getNZ(x_u);
      input(NLPSOL_LBG).getNZ(g_l);
      input(NLPSOL_UBG).getNZ(g_u);
      return true;
    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "get_bounds_info failed: " << ex.what() << endl;
      return false;
    }
  }

  bool IpoptInterface::get_starting_point(int n, bool init_x, double* x,
                                         bool init_z, double* z_L, double* z_U,
                                         int m, bool init_lambda,
                                         double* lambda) {
    try {
      bool warmstart = hasSetOption("warm_start_init_point") &&
          option("warm_start_init_point")=="yes";
      //casadi_assert_warning(init_x, "Not initializing x");
      if (warmstart) {
        //casadi_assert_warning(init_lambda, "Not initializing lambda");
        //casadi_assert_warning(init_z, "Not initializing z");
      }

      if (init_x) {
        input(NLPSOL_X0).getNZ(x);
      }

      if (init_z) {
        // Get dual solution (simple bounds)
        vector<double>& lambda_x = input(NLPSOL_LAM_X0).data();
        for (int i=0; i<lambda_x.size(); ++i) {
          z_L[i] = max(0., -lambda_x[i]);
          z_U[i] = max(0., lambda_x[i]);
        }
      }

      if (init_lambda) {
        input(NLPSOL_LAM_G0).getNZ(lambda);
      }

      return true;
    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "get_starting_point failed: " << ex.what() << endl;
      return false;
    }
  }

  void IpoptInterface::get_nlp_info(int& n, int& m, int& nnz_jac_g, int& nnz_h_lag) {
    try {
      n = nx_;               // number of variables
      m = ng_;               // number of constraints

      // Get Jacobian sparsity pattern
      if (nlp_.output(NL_G).nnz()==0)
        nnz_jac_g = 0;
      else
        nnz_jac_g = jacG().nnz_out(0);

      // Get Hessian sparsity pattern
      if (exact_hessian_)
        nnz_h_lag = hessLag().sparsity_out(0).nnz_upper();
      else
        nnz_h_lag = 0;
    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "get_nlp_info failed: " << ex.what() << endl;
    }
  }

  int IpoptInterface::get_number_of_nonlinear_variables() {
    try {
      if (!static_cast<bool>(option("pass_nonlinear_variables"))) {
        // No Hessian has been interfaced
        return -1;
      } else {
        // Number of variables that appear nonlinearily
        int nv = 0;

        // Loop over the cols
        const Sparsity& spHessLag = this->spHessLag();
        const int* colind = spHessLag.colind();
        int ncol = spHessLag.size2();
        for (int i=0; i<ncol; ++i) {
          // If the col contains any non-zeros, the corresponding variable appears nonlinearily
          if (colind[i]!=colind[i+1])
            nv++;
        }

        // Return the number
        return nv;
      }
    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "get_number_of_nonlinear_variables failed: " << ex.what() << endl;
      return -1;
    }
  }

  bool IpoptInterface::get_list_of_nonlinear_variables(int num_nonlin_vars, int* pos_nonlin_vars) {
    try {
      // Running index
      int el = 0;

      // Loop over the cols
      const Sparsity& spHessLag = this->spHessLag();
      const int* colind = spHessLag.colind();
      int ncol = spHessLag.size2();
      for (int i=0; i<ncol; ++i) {
        // If the col contains any non-zeros, the corresponding variable appears nonlinearily
        if (colind[i]!=colind[i+1]) {
          pos_nonlin_vars[el++] = i;
        }
      }

      // Assert number and return
      casadi_assert(el==num_nonlin_vars);
      return true;
    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "get_list_of_nonlinear_variables failed: " << ex.what() << endl;
      return false;
    }
  }

  bool IpoptInterface::get_var_con_metadata(int n,
                                           map<string, vector<string> >& var_string_md,
                                           map<string, vector<int> >& var_integer_md,
                                           map<string, vector<double> >& var_numeric_md,
                                           int m,
                                           map<string, vector<string> >& con_string_md,
                                           map<string, vector<int> >& con_integer_md,
                                           map<string, vector<double> >& con_numeric_md) {
    if (hasSetOption("var_string_md")) {
      Dict dict = option("var_string_md");
      for (Dict::const_iterator it=dict.begin(); it!=dict.end(); ++it) {
        string key = it->first; // Get the key
        vector<string> entry = it->second; // Get the entry
        // Check length for consistency
        casadi_assert_message(entry.size()==n, "Inconsistent length of IPOPT metadata.");
        var_string_md[key] = entry; // Save to IPOPT data structure
      }
    }

    if (hasSetOption("var_integer_md")) {
      Dict dict = option("var_integer_md");
      for (Dict::const_iterator it=dict.begin(); it!=dict.end(); ++it) {
        string key = it->first; // Get the key
        vector<int> entry = it->second; // Get the entry
        // Check length for consistency
        casadi_assert_message(entry.size()==n, "Inconsistent length of IPOPT metadata.");
        var_integer_md[key] = entry; // Save to IPOPT data structure
      }
    }

    if (hasSetOption("var_numeric_md")) {
      Dict dict = option("var_numeric_md");
      for (Dict::const_iterator it=dict.begin(); it!=dict.end(); ++it) {
        string key = it->first; // Get the key
        vector<double> entry = it->second; // Get the entry
        // Check length for consistency
        casadi_assert_message(entry.size()==n, "Inconsistent length of IPOPT metadata.");
        var_numeric_md[key] = entry; // Save to IPOPT data structure
      }
    }

    if (hasSetOption("con_string_md")) {
      Dict dict = option("con_string_md");
      for (Dict::const_iterator it=dict.begin(); it!=dict.end(); ++it) {
        string key = it->first; // Get the key
        vector<string> entry = it->second; // Get the entry
        // Check length for consistency
        casadi_assert_message(entry.size()==m, "Inconsistent length of IPOPT metadata.");
        con_string_md[key] = entry; // Save to IPOPT data structure
      }
    }

    if (hasSetOption("con_integer_md")) {
      Dict dict = option("con_integer_md");
      for (Dict::const_iterator it=dict.begin(); it!=dict.end(); ++it) {
        string key = it->first; // Get the key
        vector<int> entry = it->second; // Get the entry
        // Check length for consistency
        casadi_assert_message(entry.size()==m, "Inconsistent length of IPOPT metadata.");
        con_integer_md[key] = entry; // Save to IPOPT data structure
      }
    }

    if (hasSetOption("con_numeric_md")) {
      Dict dict = option("con_numeric_md");
      for (Dict::const_iterator it=dict.begin(); it!=dict.end(); ++it) {
        string key = it->first; // Get the key
        vector<double> entry = it->second; // Get the entry
        // Check length for consistency
        casadi_assert_message(entry.size()==m, "Inconsistent length of IPOPT metadata.");
        con_numeric_md[key] = entry; // Save to IPOPT data structure
      }
    }

    return true;
  }

  void IpoptInterface::finalize_metadata(int n,
                                        const map<string, vector<string> >& var_string_md,
                                        const map<string, vector<int> >& var_integer_md,
                                        const map<string, vector<double> >& var_numeric_md,
                                        int m,
                                        const map<string, vector<string> >& con_string_md,
                                        const map<string, vector<int> >& con_integer_md,
                                        const map<string, vector<double> >& con_numeric_md) {

    if (!var_string_md.empty()) {
      Dict dict;
      for (map<string, vector<string> >::const_iterator it=var_string_md.begin();
          it!=var_string_md.end(); ++it) {
        dict[it->first] = it->second;
      }
      stats_["var_string_md"] = dict;
    }

    if (!var_integer_md.empty()) {
      Dict dict;
      for (map<string, vector<int> >::const_iterator it=var_integer_md.begin();
          it!=var_integer_md.end(); ++it) {
        dict[it->first] = it->second;
      }
      stats_["var_integer_md"] = dict;
    }

    if (!var_numeric_md.empty()) {
      Dict dict;
      for (map<string, vector<double> >::const_iterator it=var_numeric_md.begin();
          it!=var_numeric_md.end(); ++it) {
        dict[it->first] = it->second;
      }
      stats_["var_numeric_md"] = dict;
    }

    if (!con_string_md.empty()) {
      Dict dict;
      for (map<string, vector<string> >::const_iterator it=con_string_md.begin();
          it!=con_string_md.end(); ++it) {
        dict[it->first] = it->second;
      }
      stats_["con_string_md"] = dict;
    }

    if (!con_integer_md.empty()) {
      Dict dict;
      for (map<string, vector<int> >::const_iterator it=con_integer_md.begin();
          it!=con_integer_md.end(); ++it) {
        dict[it->first] = it->second;
      }
      stats_["con_integer_md"] = dict;
    }

    if (!con_numeric_md.empty()) {
      Dict dict;
      for (map<string, vector<double> >::const_iterator it=con_numeric_md.begin();
          it!=con_numeric_md.end(); ++it) {
        dict[it->first] = it->second;
      }
      stats_["con_numeric_md"] = dict;
    }
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

  DMatrix IpoptInterface::getReducedHessian() {
#ifndef WITH_SIPOPT
    casadi_error("This feature requires sIPOPT support. Please consult the CasADi documentation.");
#else // WITH_SIPOPT
#ifndef WITH_CASADI_PATCH
    casadi_error("Retrieving the Hessian requires the CasADi sIPOPT patch. "
                 "Please consult the CasADi documentation.");
#else // WITH_CASADI_PATCH
    return red_hess_;
#endif // WITH_SIPOPT
#endif // WITH_CASADI_PATCH

  }

  // Convert a float to a string of an exact length.
  // First it tries fixed precision, then falls back to exponential notation.
  //
  // todo(jaeandersson,jgillis): needs either review or unit tests
  // because it throws exceptions if it fail.
  std::string formatFloat(double x, int totalWidth, int maxPrecision, int fallbackPrecision) {
    std::ostringstream out0;
    out0 << fixed << setw(totalWidth) << setprecision(maxPrecision) << x;
    std::string ret0 = out0.str();
    if (ret0.length() == totalWidth) {
      return ret0;
    } else if (ret0.length() > totalWidth) {
      std::ostringstream out1;
      out1 << setw(totalWidth) << setprecision(fallbackPrecision) << x;
      std::string ret1 = out1.str();
      if (ret1.length() != totalWidth)
        casadi_error(
          "ipopt timing formatting fallback is bugged, sorry about that."
          << "expected " << totalWidth <<  " digits, but got " << ret1.length()
          << ", string: \"" << ret1 << "\", number: " << x);
      return ret1;
    } else {
      casadi_error("ipopt timing formatting is bugged, sorry about that.");
    }
  }

  // Print out a beautiful timing summary.
  // Will print one row per tuple in the input vector.
  // The tuple must contain {name, # evals, total time in seconds}.
  void IpoptInterface::timingSummary(
    std::vector<std::tuple<std::string, int, DiffTime> > & xs) {
    // get padding of names
    int maxNameLen = 0;
    for (int k=0; k < xs.size(); ++k) {
      const int len = std::get<0>(xs[k]).length();
      maxNameLen = max(maxNameLen, len);
    }

    std::ostringstream out;
    std::string blankName(maxNameLen, ' ');
    out
      << blankName
      << "      user           real      num           mean             mean"
      << endl << blankName
      << "      time           time     evals       user time        real time"
      << endl;
    for (int k=0; k < xs.size(); ++k) {
      const std::tuple<std::string, int, DiffTime> x = xs[k];
      const std::string name = std::get<0>(x);
      const int n = std::get<1>(x);
      DiffTime dt = std::get<2>(x);

      // don't print out this row if there were 0 calls
      if (n == 0)
        continue;

      out
        << setw(maxNameLen) << name << " "
        << formatFloat(dt.user, 9, 3, 3) << " [s]  "
        << formatFloat(dt.real, 9, 3, 3) << " [s]";
      if (n == -1) {
        // things like main loop don't have # evals
        out << endl;
      } else {
        out
          << " "
          << setw(5) << n;
        if (n < 2) {
          out << endl;
        } else {
          // only print averages if there is more than 1 eval
          out
            << " "
            << formatFloat(1000.0*dt.user/n, 10, 2, 3) << " [ms]  "
            << formatFloat(1000.0*dt.real/n, 10, 2, 3) << " [ms]"
            << endl;
        }
      }
    }
    // I'm worried that the following will set stdout stream state:
    // userOut() << out;

    // Just convert it to a string first to be safe.
    // todo(jaeandersson/jgillis): please review.
    const std::string ret = out.str();
    userOut() << ret;
  }

  void IpoptInterface::set_x(const double *x) {
    // Is a recalculation needed
    if (new_x_ || !equal(x, x+nx_, xk_)) {
      copy_n(x, nx_, xk_);
      have_fk_ = have_gk_ = have_hess_lk_ = have_grad_lk_ = have_grad_fk_
        = have_jac_gk_ = false;
      new_x_ = false;
    }
  }

  void IpoptInterface::set_lam_f(double lam_f) {
    // Is a recalculation needed
    if (new_lam_f_ || lam_f != lam_fk_) {
      lam_fk_ = lam_f;
      have_hess_lk_ = have_grad_lk_ = false;
      new_lam_f_ = false;
    }
  }

  void IpoptInterface::set_lam_g(const double *lam_g) {
    // Is a recalculation needed
    if (new_lam_g_ || !equal(lam_g, lam_g+ng_, lam_gk_)) {
      copy_n(lam_g, ng_, lam_gk_);
      have_hess_lk_ = have_grad_lk_ = false;
      new_lam_g_ = false;
    }
  }

  int IpoptInterface::calc_f(double* f) {
    // Respond to a possible Crl+C signals
    InterruptHandler::check();

    // Calculate, if needed
    if (!have_fk_) {
      fill_n(arg_, f_fcn_.n_in(), nullptr);
      arg_[F_X] = xk_;
      arg_[F_P] = p_;
      fill_n(res_, f_fcn_.n_out(), nullptr);
      res_[F_F] = &fk_;
      n_calc_f_ += 1;
      auto t_start = chrono::system_clock::now(); // start timer
      try {
        f_fcn_(arg_, res_, iw_, w_, 0);
      } catch(exception& ex) {
        // Fatal error
        userOut<true, PL_WARN>() << name() << ":calc_f failed:" << ex.what() << endl;
        return 1;
      }
      auto t_stop = chrono::system_clock::now(); // stop timer

      // Make sure not NaN or Inf
      if (!isfinite(fk_)) {
        userOut<true, PL_WARN>() << name() << ":calc_f failed: Inf or NaN detected" << endl;
        return -1;
      }

      // Update stats
      n_calc_f_ += 1;
      t_calc_f_ += chrono::duration<double>(t_stop - t_start).count();
      have_fk_ = true;
    }

    // Return to user
    if (f) *f = fk_;

    // Success
    return 0;
  }

  template<typename M>
  void IpoptInterface::setup_f() {
    const Problem<M>& nlp = nlp2_;
    vector<M> arg(F_NUM_IN);
    arg[F_X] = nlp.in[NL_X];
    arg[F_P] = nlp.in[NL_P];
    vector<M> res(F_NUM_OUT);
    res[F_F] = nlp.out[NL_F];
    f_fcn_ = Function("nlp_f", arg, res);
    alloc(f_fcn_);
  }

  int IpoptInterface::calc_g(double* g) {
    // Respond to a possible Crl+C signals
    InterruptHandler::check();

    // Calculate, if needed
    if (!have_gk_) {
      // Evaluate User function
      fill_n(arg_, g_fcn_.n_in(), nullptr);
      arg_[G_X] = xk_;
      arg_[G_P] = p_;
      fill_n(res_, g_fcn_.n_out(), nullptr);
      res_[G_G] = gk_;
      auto t_start = chrono::system_clock::now(); // start timer
      try {
        g_fcn_(arg_, res_, iw_, w_, 0);
      } catch(exception& ex) {
        // Fatal error
        userOut<true, PL_WARN>() << name() << ":calc_g failed:" << ex.what() << endl;
        return 1;
      }
      auto t_stop = chrono::system_clock::now(); // stop timer

      // Make sure not NaN or Inf
      if (!all_of(gk_, gk_+ng_, [](double v) { return isfinite(v);})) {
        userOut<true, PL_WARN>() << name() << ":calc_g failed: NaN or Inf detected" << endl;
        return -1;
      }

      // Update stats
      n_calc_g_ += 1;
      t_calc_g_ += chrono::duration<double>(t_stop - t_start).count();
      have_gk_ = true;
    }

    // Return to user
    casadi_copy(gk_, ng_, g);

    // Success
    return 0;
  }

  template<typename M>
  void IpoptInterface::setup_g() {
    const Problem<M>& nlp = nlp2_;
    vector<M> arg(G_NUM_IN);
    arg[G_X] = nlp.in[NL_X];
    arg[G_P] = nlp.in[NL_P];
    vector<M> res(G_NUM_OUT);
    res[G_G] = nlp.out[NL_G];
    g_fcn_ = Function("nlp_g", arg, res);
    alloc(g_fcn_);
  }

  int IpoptInterface::calc_grad_f(double* grad_f) {
    // Respond to a possible Crl+C signals
    InterruptHandler::check();

    // Calculate, if needed
    if (!have_grad_fk_) {
      // Evaluate User function
      fill_n(arg_, grad_f_fcn_.n_in(), nullptr);
      arg_[GF_X] = xk_;
      arg_[GF_P] = p_;
      fill_n(res_, grad_f_fcn_.n_out(), nullptr);
      res_[GF_GF] = grad_fk_;
      auto t_start = chrono::system_clock::now(); // start timer
      try {
        grad_f_fcn_(arg_, res_, iw_, w_, 0);
      } catch(exception& ex) {
        // Fatal error
        userOut<true, PL_WARN>() << name() << ":calc_grad_f failed:" << ex.what() << endl;
        return 1;
      }
      auto t_stop = chrono::system_clock::now(); // stop timer

      // Make sure not NaN or Inf
      if (!all_of(grad_fk_, grad_fk_+nx_, [](double v) { return isfinite(v);})) {
        userOut<true, PL_WARN>() << name() << ":calc_grad_f failed: NaN or Inf detected" << endl;
        return -1;
      }

      // Update stats
      n_calc_grad_f_ += 1;
      t_calc_grad_f_ += chrono::duration<double>(t_stop - t_start).count();
      have_grad_fk_ = true;
    }

    // Return to user
    casadi_copy(grad_fk_, nx_, grad_f);

    // Success
    return 0;
  }

  template<typename M>
  void IpoptInterface::setup_grad_f() {
    const Problem<M>& nlp = nlp2_;
    vector<M> arg(GF_NUM_IN);
    arg[GF_X] = nlp.in[NL_X];
    arg[GF_P] = nlp.in[NL_P];
    vector<M> res(GF_NUM_OUT);
    res[GF_GF] = M::gradient(nlp.out[NL_F], nlp.in[NL_X]);
    grad_f_fcn_ = Function("nlp_grad_f", arg, res);
    alloc(grad_f_fcn_);
  }

  int IpoptInterface::calc_jac_g(double* jac_g) {
    // Respond to a possible Crl+C signals
    InterruptHandler::check();

    // Calculate, if needed
    if (!have_jac_gk_) {
      // Evaluate User function
      fill_n(arg_, jac_g_fcn_.n_in(), nullptr);
      arg_[JG_X] = xk_;
      arg_[JG_P] = p_;
      fill_n(res_, jac_g_fcn_.n_out(), nullptr);
      res_[JG_JG] = jac_gk_;
      auto t_start = chrono::system_clock::now(); // start timer
      try {
        jac_g_fcn_(arg_, res_, iw_, w_, 0);
      } catch(exception& ex) {
        // Fatal error
        userOut<true, PL_WARN>() << name() << ":calc_jac_g failed:" << ex.what() << endl;
        return 1;
      }
      auto t_stop = chrono::system_clock::now(); // stop timer

      // Make sure not NaN or Inf
      if (!all_of(jac_gk_, jac_gk_+jacg_sp_.nnz(), [](double v) { return isfinite(v);})) {
        userOut<true, PL_WARN>() << name() << ":calc_jac_g failed: NaN or Inf detected" << endl;
        return -1;
      }

      // Update stats
      n_calc_jac_g_ += 1;
      t_calc_jac_g_ += chrono::duration<double>(t_stop - t_start).count();
      have_jac_gk_ = true;
    }

    // Return to user
    casadi_copy(jac_gk_, jacg_sp_.nnz(), jac_g);

    // Success
    return 0;
  }

  template<typename M>
  void IpoptInterface::setup_jac_g() {
    const Problem<M>& nlp = nlp2_;
    vector<M> arg(JG_NUM_IN);
    arg[JG_X] = nlp.in[NL_X];
    arg[JG_P] = nlp.in[NL_P];
    vector<M> res(JG_NUM_OUT);
    res[JG_JG] = M::jacobian(nlp.out[NL_G], nlp.in[NL_X]);
    jac_g_fcn_ = Function("nlp_jac_g", arg, res);
    jacg_sp_ = res[JG_JG].sparsity();
    alloc(jac_g_fcn_);
  }

  int IpoptInterface::calc_hess_l(double* hess_l) {
    // Respond to a possible Crl+C signals
    InterruptHandler::check();

    // Calculate, if needed
    if (!have_hess_lk_) {
      // Evaluate User function
      fill_n(arg_, hess_l_fcn_.n_in(), nullptr);
      arg_[HL_X] = xk_;
      arg_[HL_P] = p_;
      arg_[HL_LAM_F] = &lam_fk_;
      arg_[HL_LAM_G] = lam_gk_;
      fill_n(res_, hess_l_fcn_.n_out(), nullptr);
      res_[HL_HL] = hess_lk_;
      auto t_start = chrono::system_clock::now(); // start timer
      try {
        hess_l_fcn_(arg_, res_, iw_, w_, 0);
      } catch(exception& ex) {
        // Fatal error
        userOut<true, PL_WARN>() << name() << ":calc_hess_l failed:" << ex.what() << endl;
        return 1;
      }
      auto t_stop = chrono::system_clock::now(); // stop timer

      // Make sure not NaN or Inf
      if (!all_of(hess_lk_, hess_lk_+hesslag_sp_.nnz(), [](double v) { return isfinite(v);})) {
        userOut<true, PL_WARN>() << name() << ":calc_hess_l failed: NaN or Inf detected" << endl;
        return -1;
      }

      // Update stats
      n_calc_hess_l_ += 1;
      t_calc_hess_l_ += chrono::duration<double>(t_stop - t_start).count();
      have_hess_lk_ = true;
    }

    // Return to user
    casadi_copy(hess_lk_, hesslag_sp_.nnz(), hess_l);

    // Success
    return 0;
  }

  template<typename M>
  void IpoptInterface::setup_hess_l() {
    const Problem<M>& nlp = nlp2_;
    vector<M> arg(HL_NUM_IN);
    M x = arg[HL_X] = nlp.in[NL_X];
    arg[HL_P] = nlp.in[NL_P];
    M f = nlp.out[NL_F];
    M g = nlp.out[NL_G];
    M lam_f = arg[HL_LAM_F] = M::sym("lam_f", f.sparsity());
    M lam_g = arg[HL_LAM_G] = M::sym("lam_g", g.sparsity());
    vector<M> res(HL_NUM_OUT);
    res[HL_HL] = M::hessian(dot(lam_f, f) + dot(lam_g, g), x);
    hess_l_fcn_ = Function("nlp_hess_l", arg, res);
    hesslag_sp_ = res[HL_HL].sparsity();
    alloc(hess_l_fcn_);
  }

  template<typename M>
  void IpoptInterface::setup() {
    // Objective
    setup_f<M>();

    // Constraints
    setup_g<M>();

    // Gradient of the objective
    setup_grad_f<M>();

    // Jacobian of the constraits
    setup_jac_g<M>();

    // Hessian of the Lagrangian
    if (exact_hessian_) {
      setup_hess_l<M>();
    }
  }

} // namespace casadi
