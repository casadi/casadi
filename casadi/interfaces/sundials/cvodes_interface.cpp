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


#include "cvodes_interface.hpp"

#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/function/linear_solver_internal.hpp"

using namespace std;

namespace casadi {

  extern "C"
  int CASADI_INTEGRATOR_CVODES_EXPORT
      casadi_register_integrator_cvodes(Integrator::Plugin* plugin) {
    plugin->creator = CvodesInterface::creator;
    plugin->name = "cvodes";
    plugin->doc = CvodesInterface::meta_doc.c_str();;
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_INTEGRATOR_CVODES_EXPORT casadi_load_integrator_cvodes() {
    Integrator::registerPlugin(casadi_register_integrator_cvodes);
  }

  CvodesInterface::CvodesInterface(const std::string& name, const XProblem& dae)
    : SundialsInterface(name, dae) {

    addOption("linear_multistep_method",          OT_STRING,              "bdf",
              "Integrator scheme", "bdf|adams");
    addOption("nonlinear_solver_iteration",       OT_STRING,              "newton",
              "", "newton|functional");
    addOption("fsens_all_at_once",                OT_BOOLEAN,             true,
              "Calculate all right hand sides of the sensitivity equations at once");
    addOption("disable_internal_warnings",        OT_BOOLEAN,             false,
              "Disable CVodes internal warning messages");
    addOption("monitor",                          OT_STRINGVECTOR,        GenericType(),
              "", "res|resB|resQB|reset|psetupB|djacB", true);

    mem_ = 0;

    x0_ = x_ = q_ = 0;
    rx0_ = rx_ = rq_ = 0;

    isInitAdj_ = false;
    disable_internal_warnings_ = false;
  }

  void CvodesInterface::freeCVodes() {
    if (mem_) { CVodeFree(&mem_); mem_ = 0;}

    // Forward integration
    if (x0_) { N_VDestroy_Serial(x0_); x0_ = 0; }
    if (x_) { N_VDestroy_Serial(x_); x_ = 0; }
    if (q_) { N_VDestroy_Serial(q_); q_ = 0; }

    // Backward integration
    if (rx0_) { N_VDestroy_Serial(rx0_); rx0_ = 0; }
    if (rx_)  { N_VDestroy_Serial(rx_);  rx_  = 0; }
    if (rq_)  { N_VDestroy_Serial(rq_);  rq_  = 0; }

    // Sensitivities of the forward integration
    for (vector<N_Vector>::iterator it=xF0_.begin(); it != xF0_.end(); ++it)
        if (*it) { N_VDestroy_Serial(*it); *it=0;}
    for (vector<N_Vector>::iterator it=xF_.begin(); it != xF_.end(); ++it)
        if (*it) { N_VDestroy_Serial(*it); *it=0;}
    for (vector<N_Vector>::iterator it=qF_.begin(); it != qF_.end(); ++it)
        if (*it) { N_VDestroy_Serial(*it); *it=0;}
  }

  CvodesInterface::~CvodesInterface() {
    freeCVodes();
  }

  void CvodesInterface::init() {
    log("CvodesInterface::init", "begin");

    // Free memory if already initialized
    freeCVodes();

    // Initialize the base classes
    SundialsInterface::init();

    // Reset checkpoints counter
    ncheck_ = 0;

    // Algebraic variables not supported
    casadi_assert_message(nz_==0 && nrz_==0,
                          "CVODES does not support algebraic variables");

    // Read options
    monitor_rhsB_  = monitored("resB");
    monitor_rhs_   = monitored("res");
    monitor_rhsQB_ = monitored("resQB");

    // Sundials return flag
    int flag;

    if (option("linear_multistep_method")=="adams")  lmm_ = CV_ADAMS;
    else if (option("linear_multistep_method")=="bdf") lmm_ = CV_BDF;
    else
      throw CasadiException("Unknown linear multistep method");

    if (option("nonlinear_solver_iteration")=="newton") iter_ = CV_NEWTON;
    else if (option("nonlinear_solver_iteration")=="functional") iter_ = CV_FUNCTIONAL;
    else
      throw CasadiException("Unknown nonlinear solver iteration");

    // Create CVodes memory block
    mem_ = CVodeCreate(lmm_, iter_);
    if (mem_==0) throw CasadiException("CVodeCreate: Creation failed");

    // Allocate n-vectors for ivp
    x0_ = N_VMake_Serial(nx_, x0().ptr());
    x_ = N_VMake_Serial(nx_, xf().ptr());

    // Disable internal warning messages?
    disable_internal_warnings_ = option("disable_internal_warnings");

    // Set error handler function
    flag = CVodeSetErrHandlerFn(mem_, ehfun_wrapper, this);
    if (flag != CV_SUCCESS) cvodes_error("CVodeSetErrHandlerFn", flag);

    // Initialize CVodes
    double t0 = 0;
    flag = CVodeInit(mem_, rhs_wrapper, t0, x0_);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeInit", flag);

    // Set tolerances
    flag = CVodeSStolerances(mem_, reltol_, abstol_);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeInit", flag);

    // Maximum number of steps
    CVodeSetMaxNumSteps(mem_, option("max_num_steps").toInt());
    if (flag != CV_SUCCESS) cvodes_error("CVodeSetMaxNumSteps", flag);

    // attach a linear solver
    switch (linsol_f_) {
    case SD_DENSE:
      initDenseLinearSolver();
      break;
    case SD_BANDED:
      initBandedLinearSolver();
      break;
    case SD_ITERATIVE:
      initIterativeLinearSolver();
      break;
    case SD_USER_DEFINED:
      initUserDefinedLinearSolver();
      break;
    }

    // Set user data
    flag = CVodeSetUserData(mem_, this);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeSetUserData", flag);

    // Quadrature equations
    if (nq_>0) {
      // Allocate n-vectors for quadratures
      q_ = N_VMake_Serial(nq_, qf().ptr());

      // Initialize quadratures in CVodes
      N_VConst(0.0, q_);
      flag = CVodeQuadInit(mem_, rhsQ_wrapper, q_);
      if (flag != CV_SUCCESS) cvodes_error("CVodeQuadInit", flag);

      // Should the quadrature errors be used for step size control?
      if (option("quad_err_con").toInt()) {
        flag = CVodeSetQuadErrCon(mem_, true);
        if (flag != CV_SUCCESS) cvodes_error("CVodeSetQuadErrCon", flag);

        // Quadrature error tolerances
        // TODO(Joel): vector absolute tolerances
        flag = CVodeQuadSStolerances(mem_, reltol_, abstol_);
        if (flag != CV_SUCCESS) cvodes_error("CVodeQuadSStolerances", flag);
      }
    }

    // Adjoint sensitivity problem
    if (!g_.isNull()) {

      // Allocate n-vectors for backward integration
      rx0_ = N_VMake_Serial(nrx_, rx0().ptr());
      rx_ = N_VMake_Serial(nrx_, rxf().ptr());
      rq_ = N_VMake_Serial(nrq_, rqf().ptr());

      // Get the number of steos per checkpoint
      int Nd = option("steps_per_checkpoint");

      // Get the interpolation type
      int interpType;
      if (option("interpolation_type")=="hermite")
        interpType = CV_HERMITE;
      else if (option("interpolation_type")=="polynomial")
        interpType = CV_POLYNOMIAL;
      else
        throw CasadiException("\"interpolation_type\" must be \"hermite\" or \"polynomial\"");

      // Initialize adjoint sensitivities
      flag = CVodeAdjInit(mem_, Nd, interpType);
      if (flag != CV_SUCCESS) cvodes_error("CVodeAdjInit", flag);

      isInitAdj_ = false;
    }
  }


  void CvodesInterface::initAdj() {

    // Create backward problem (use the same lmm and iter)
    int flag = CVodeCreateB(mem_, lmm_, iter_, &whichB_);
    if (flag != CV_SUCCESS) cvodes_error("CVodeCreateB", flag);

    // Initialize the backward problem
    double tB0 = grid_.back();
    flag = CVodeInitB(mem_, whichB_, rhsB_wrapper, tB0, rx0_);
    if (flag != CV_SUCCESS) cvodes_error("CVodeInitB", flag);

    // Set tolerances
    flag = CVodeSStolerancesB(mem_, whichB_, reltolB_, abstolB_);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeSStolerancesB", flag);

    // User data
    flag = CVodeSetUserDataB(mem_, whichB_, this);
    if (flag != CV_SUCCESS) cvodes_error("CVodeSetUserDataB", flag);

    // attach linear solver to backward problem
    switch (linsol_g_) {
    case SD_DENSE:
      initDenseLinearSolverB();
      break;
    case SD_BANDED:
      initBandedLinearSolverB();
      break;
    case SD_ITERATIVE:
      initIterativeLinearSolverB();
      break;
    case SD_USER_DEFINED:
      initUserDefinedLinearSolverB();
      break;
    }

    // Quadratures for the backward problem
    N_VConst(0.0, rq_);
    flag = CVodeQuadInitB(mem_, whichB_, rhsQB_wrapper, rq_);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeQuadInitB", flag);

    if (option("quad_err_con").toInt()) {
      flag = CVodeSetQuadErrConB(mem_, whichB_, true);
      if (flag != CV_SUCCESS) cvodes_error("CVodeSetQuadErrConB", flag);

      flag = CVodeQuadSStolerancesB(mem_, whichB_, reltolB_, abstolB_);
      if (flag != CV_SUCCESS) cvodes_error("CVodeQuadSStolerancesB", flag);
    }

    // Mark initialized
    isInitAdj_ = true;
  }

  void CvodesInterface::rhs(double t, const double* x, double* xdot) {
    log("CvodesInterface::rhs", "begin");

    // Get time
    time1 = clock();

    // Pass input
    f_.setInputNZ(&t, DAE_T);
    f_.setInputNZ(x, DAE_X);
    f_.setInput(p(), DAE_P);

    if (monitor_rhs_) {
      userOut() << "t       = " << t << endl;
      userOut() << "x       = " << f_.input(DAE_X) << endl;
      userOut() << "p       = " << f_.input(DAE_P) << endl;
    }
    // Evaluate
    f_.evaluate();

    if (monitor_rhs_) {
      userOut() << "xdot       = " << f_.output(DAE_ODE)<< endl;
    }

    // Get results
    f_.getOutputNZ(xdot);

    // Log time
    time2 = clock();
    t_res += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;

    log("CvodesInterface::rhs", "end");

  }

  int CvodesInterface::rhs_wrapper(double t, N_Vector x, N_Vector xdot, void *user_data) {
    try {
      casadi_assert(user_data);
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      this_->rhs(t, NV_DATA_S(x), NV_DATA_S(xdot));
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhs failed: " << e.what() << endl;
      return 1;
    }
  }

  void CvodesInterface::reset() {
    casadi_msg("CvodesInterface::reset begin");

    // Reset the base classes
    SundialsInterface::reset();

    if (monitored("reset")) {
      userOut() << "initial state: " << endl;
      userOut() << "p = " << p() << endl;
      userOut() << "x0 = " << x0() << endl;
    }

    // Reset timers
    t_res = t_fres = t_jac = t_lsolve = t_lsetup_jac = t_lsetup_fac = 0;

    // Re-initialize
    int flag = CVodeReInit(mem_, grid_.front(), x0_);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeReInit", flag);

    // Re-initialize quadratures
    if (nq_>0) {
      N_VConst(0.0, q_);
      flag = CVodeQuadReInit(mem_, q_);
      if (flag != CV_SUCCESS) cvodes_error("CVodeQuadReInit", flag);
    }

    // Turn off sensitivities
    flag = CVodeSensToggleOff(mem_);
    if (flag != CV_SUCCESS) cvodes_error("CVodeSensToggleOff", flag);

    // Re-initialize backward integration
    if (nrx_>0) {
      flag = CVodeAdjReInit(mem_);
      if (flag != CV_SUCCESS) cvodes_error("CVodeAdjReInit", flag);
    }

    // Set the stop time of the integration -- don't integrate past this point
    if (stop_at_end_) setStopTime(grid_.back());
    casadi_msg("CvodesInterface::reset end");
  }

  void CvodesInterface::advance(int k) {
    double t_out = grid_.at(k);

    casadi_msg("CvodesInterface::integrate(" << t_out << ") begin");

    casadi_assert_message(t_out>=grid_.front(),
                          "CvodesInterface::integrate(" << t_out << "): "
                          "Cannot integrate to a time earlier than t0 ("
                          << grid_.front() << ")");
    casadi_assert_message(t_out<=grid_.back() || !stop_at_end_, "CvodesInterface::integrate("
                          << t_out << "):"
                          " Cannot integrate past a time later than tf (" << grid_.back() << ") "
                          "unless stop_at_end is set to False.");

    int flag;

    // tolerance
    double ttol = 1e-9;
    if (fabs(t_-t_out)<ttol) {
      return;
    }
    if (nrx_>0) {
      flag = CVodeF(mem_, t_out, x_, &t_, CV_NORMAL, &ncheck_);
      if (flag!=CV_SUCCESS && flag!=CV_TSTOP_RETURN) cvodes_error("CVodeF", flag);

    } else {
      flag = CVode(mem_, t_out, x_, &t_, CV_NORMAL);
      if (flag!=CV_SUCCESS && flag!=CV_TSTOP_RETURN) cvodes_error("CVode", flag);
    }

    if (nq_>0) {
      double tret;
      flag = CVodeGetQuad(mem_, &tret, q_);
      if (flag!=CV_SUCCESS) cvodes_error("CVodeGetQuad", flag);
    }

    // Print statistics
    if (option("print_stats")) printStats(userOut());

    if (gather_stats_) {
      long nsteps, nfevals, nlinsetups, netfails;
      int qlast, qcur;
      double hinused, hlast, hcur, tcur;
      int flag = CVodeGetIntegratorStats(mem_, &nsteps, &nfevals, &nlinsetups, &netfails, &qlast,
                                         &qcur, &hinused, &hlast, &hcur, &tcur);
      if (flag!=CV_SUCCESS) cvodes_error("CVodeGetIntegratorStats", flag);

      stats_["nsteps"] = 1.0*nsteps;
      stats_["nlinsetups"] = 1.0*nlinsetups;

    }

    casadi_msg("CvodesInterface::integrate(" << t_out << ") end");
  }

  void CvodesInterface::resetB() {
    casadi_msg("CvodesInterface::resetB begin");

    // Reset the base classes
    SundialsInterface::resetB();

    int flag;
    if (isInitAdj_) {

      flag = CVodeReInitB(mem_, whichB_, grid_.back(), rx0_);
      if (flag != CV_SUCCESS) cvodes_error("CVodeReInitB", flag);

      N_VConst(0.0, rq_);
      flag = CVodeQuadReInitB(mem_, whichB_, rq_);
      if (flag!=CV_SUCCESS) cvodes_error("CVodeQuadReInitB", flag);

    } else {
      // Initialize the adjoint integration
      initAdj();
    }
    casadi_msg("CvodesInterface::resetB end");
  }

  void CvodesInterface::retreat(int k) {
    double t_out = grid_.at(k);

    casadi_msg("CvodesInterface::retreat(" << t_out << ") begin");
    int flag;

    // Integrate backward to t_out
    flag = CVodeB(mem_, t_out, CV_NORMAL);
    if (flag<CV_SUCCESS) cvodes_error("CVodeB", flag);

    // Get the sensitivities
    double tret;
    flag = CVodeGetB(mem_, whichB_, &tret, rx_);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeGetB", flag);

    flag = CVodeGetQuadB(mem_, whichB_, &tret, rq_);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeGetQuadB", flag);


    if (gather_stats_) {
      long nsteps, nfevals, nlinsetups, netfails;
      int qlast, qcur;
      double hinused, hlast, hcur, tcur;
      CVodeMem cv_mem = static_cast<CVodeMem>(mem_);
      CVadjMem ca_mem = cv_mem->cv_adj_mem;
      CVodeBMem cvB_mem = ca_mem->cvB_mem;

      int flag = CVodeGetIntegratorStats(cvB_mem->cv_mem, &nsteps,
                                         &nfevals, &nlinsetups, &netfails, &qlast, &qcur,
                                         &hinused, &hlast, &hcur, &tcur);
      if (flag!=CV_SUCCESS) cvodes_error("CVodeGetIntegratorStats", flag);

      stats_["nstepsB"] = 1.0*nsteps;
      stats_["nlinsetupsB"] = 1.0*nlinsetups;

    }
    casadi_msg("CvodesInterface::retreat(" << t_out << ") end");
  }

  void CvodesInterface::printStats(std::ostream &stream) const {
    long nsteps, nfevals, nlinsetups, netfails;
    int qlast, qcur;
    double hinused, hlast, hcur, tcur;
    int flag = CVodeGetIntegratorStats(mem_, &nsteps, &nfevals, &nlinsetups, &netfails, &qlast,
                                       &qcur, &hinused, &hlast, &hcur, &tcur);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeGetIntegratorStats", flag);

    // Get the number of right hand side evaluations in the linear solver
    long nfevals_linsol=0;
    switch (linsol_f_) {
    case SD_DENSE:
    case SD_BANDED:
      flag = CVDlsGetNumRhsEvals(mem_, &nfevals_linsol);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsGetNumRhsEvals", flag);
      break;
    case SD_ITERATIVE:
      flag = CVSpilsGetNumRhsEvals(mem_, &nfevals_linsol);
      if (flag!=CV_SUCCESS) cvodes_error("CVSpilsGetNumRhsEvals", flag);
      break;
    default:
      nfevals_linsol = 0;
    }

    stream << "number of steps taken by CVODES:          " << nsteps << std::endl;
    stream << "number of calls to the user's f function: " << (nfevals + nfevals_linsol)
           << std::endl;
    stream << "   step calculation:                      " << nfevals << std::endl;
    stream << "   linear solver:                         " << nfevals_linsol << std::endl;
    stream << "number of calls made to the linear solver setup function: " << nlinsetups
           << std::endl;
    stream << "number of error test failures: " << netfails << std::endl;
    stream << "method order used on the last internal step: " << qlast << std::endl;
    stream << "method order to be used on the next internal step: " << qcur << std::endl;
    stream << "actual value of initial step size: " << hinused << std::endl;
    stream << "step size taken on the last internal step: " << hlast << std::endl;
    stream << "step size to be attempted on the next internal step: " << hcur << std::endl;
    stream << "current internal time reached: " << tcur << std::endl;
    stream << std::endl;

    stream << "number of checkpoints stored: " << ncheck_ << endl;
    stream << std::endl;

    stream << "Time spent in the ODE residual: " << t_res << " s." << endl;
    stream << "Time spent in the forward sensitivity residual: " << t_fres << " s." << endl;
    stream << "Time spent in the jacobian function or jacobian times vector function: "
           << t_jac << " s." << endl;
    stream << "Time spent in the linear solver solve function: " << t_lsolve << " s." << endl;
    stream << "Time spent to generate the jacobian in the linear solver setup function: "
           << t_lsetup_jac << " s." << endl;
    stream << "Time spent to factorize the jacobian in the linear solver setup function: "
           << t_lsetup_fac << " s." << endl;
    stream << std::endl;
  }

  map<int, string> CvodesInterface::calc_flagmap() {
    map<int, string> f;
    f[CV_SUCCESS] = "CV_SUCCESS";
    f[CV_TSTOP_RETURN] = "CV_TSTOP_RETURN";
    f[CV_ROOT_RETURN] = "CV_ROOT_RETURN";
    f[CV_WARNING] = "CV_WARNING";
    f[CV_WARNING] = "CV_WARNING";
    f[CV_TOO_MUCH_WORK] = "CV_TOO_MUCH_WORK";
    f[CV_TOO_MUCH_ACC] = "CV_TOO_MUCH_ACC";
    f[CV_ERR_FAILURE] = "CV_ERR_FAILURE";
    f[CV_CONV_FAILURE] = "CV_CONV_FAILURE";
    f[CV_LINIT_FAIL] = "CV_LINIT_FAIL";
    f[CV_LSETUP_FAIL] = "CV_LSETUP_FAIL";
    f[CV_LSOLVE_FAIL] = "CV_LSOLVE_FAIL";
    f[CV_RHSFUNC_FAIL] = "CV_RHSFUNC_FAIL";
    f[CV_FIRST_RHSFUNC_ERR] = "CV_FIRST_RHSFUNC_ERR";
    f[CV_UNREC_RHSFUNC_ERR] = "CV_UNREC_RHSFUNC_ERR";
    f[CV_RTFUNC_FAIL] = "CV_RTFUNC_FAIL";
    f[CV_MEM_FAIL] = "CV_MEM_FAIL";
    f[CV_ILL_INPUT] = "CV_ILL_INPUT";
    f[CV_NO_MALLOC] = "CV_NO_MALLOC";
    f[CV_BAD_K] = "CV_BAD_K";
    f[CV_BAD_T] = "CV_BAD_T";
    f[CV_BAD_DKY] = "CV_BAD_DKY";
    f[CV_TOO_CLOSE] = "CV_TOO_CLOSE";
    f[CV_QRHSFUNC_FAIL] = "CV_QRHSFUNC_FAIL";
    f[CV_FIRST_QRHSFUNC_ERR] = "CV_FIRST_QRHSFUNC_ERR";
    f[CV_REPTD_QRHSFUNC_ERR] = "CV_REPTD_QRHSFUNC_ERR";
    f[CV_UNREC_QRHSFUNC_ERR] = "CV_UNREC_QRHSFUNC_ERR";
    f[CV_NO_SENS ] = "CV_NO_SENS ";
    f[CV_SRHSFUNC_FAIL] = "CV_SRHSFUNC_FAIL";
    return f;
  }

  map<int, string> CvodesInterface::flagmap = CvodesInterface::calc_flagmap();

  void CvodesInterface::cvodes_error(const string& module, int flag) {
    // Find the error
    map<int, string>::const_iterator it = flagmap.find(flag);

    stringstream ss;
    if (it == flagmap.end()) {
      ss << "Unknown error (" << flag << ") from module \"" << module << "\".";
    } else {
      ss << "Module \"" << module << "\" returned flag \"" << it->second << "\".";
    }
    ss << " Consult Cvodes documentation.";
    casadi_error(ss.str());
  }

  void CvodesInterface::ehfun_wrapper(int error_code, const char *module, const char *function,
                                     char *msg, void *user_data) {
    try {
      casadi_assert(user_data);
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      this_->ehfun(error_code, module, function, msg);
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "ehfun failed: " << e.what() << endl;
    }
  }

  void CvodesInterface::ehfun(int error_code, const char *module, const char *function, char *msg) {
    if (!disable_internal_warnings_) {
      userOut<true, PL_WARN>() << msg << endl;
    }
  }

  void CvodesInterface::rhsS(int Ns, double t, N_Vector x, N_Vector xdot, N_Vector *xF,
                            N_Vector *xdotF, N_Vector tmp1, N_Vector tmp2) {
    //    casadi_assert(Ns==nfdir_);

    // Record the current cpu time
    time1 = clock();

    // Commented out since a new implementation currently cannot be tested
    casadi_error("Commented out, #884, #794.");

    // Record timings
    time2 = clock();
    t_fres += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
  }

  int CvodesInterface::rhsS_wrapper(int Ns, double t, N_Vector x, N_Vector xdot, N_Vector *xF,
                                   N_Vector *xdotF, void *user_data, N_Vector tmp1, N_Vector tmp2) {
    try {
      casadi_assert(user_data);
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      this_->rhsS(Ns, t, x, xdot, xF, xdotF, tmp1, tmp2);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "fs failed: " << e.what() << endl;
      return 1;
    }
  }

  void CvodesInterface::rhsS1(int Ns, double t, N_Vector x, N_Vector xdot, int iS, N_Vector xF,
                             N_Vector xdotF, N_Vector tmp1, N_Vector tmp2) {

    // Commented out since a new implementation currently cannot be tested
    casadi_error("Commented out, #884, #794.");
  }

  int CvodesInterface::rhsS1_wrapper(int Ns, double t, N_Vector x, N_Vector xdot, int iS,
                                    N_Vector xF, N_Vector xdotF, void *user_data,
                                    N_Vector tmp1, N_Vector tmp2) {
    try {
      casadi_assert(user_data);
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      this_->rhsS1(Ns, t, x, xdot, iS, xF, xdotF, tmp1, tmp2);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "fs failed: " << e.what() << endl;
      return 1;
    }
  }

  int CvodesInterface::rhsQ_wrapper(double t, N_Vector x, N_Vector qdot, void *user_data) {
    try {
      casadi_assert(user_data);
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      this_->rhsQ(t, NV_DATA_S(x), NV_DATA_S(qdot));
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhsQ failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::rhsQ(double t, const double* x, double* qdot) {
    // Pass input
    f_.setInputNZ(&t, DAE_T);
    f_.setInputNZ(x, DAE_X);
    f_.setInput(p(), DAE_P);

    // Evaluate
    f_.evaluate();

    // Get results
    f_.getOutputNZ(qdot, DAE_QUAD);
  }

  void CvodesInterface::rhsQS(int Ns, double t, N_Vector x, N_Vector *xF, N_Vector qdot,
                             N_Vector *qdotF, N_Vector tmp1, N_Vector tmp2) {
    //    casadi_assert(Ns==nfdir_);

    // Commented out since a new implementation currently cannot be tested
    casadi_error("Commented out, #884, #794.");
  }

  int CvodesInterface::rhsQS_wrapper(int Ns, double t, N_Vector x, N_Vector *xF, N_Vector qdot,
                                    N_Vector *qdotF, void *user_data,
                                    N_Vector tmp1, N_Vector tmp2) {
    try {
      //    casadi_assert(user_data);
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      if (!this_) {
        // SUNDIALS BUG!!!
        for (int i=0; i<Ns; ++i) N_VConst(0.0, qdotF[i]);
        return 0;
      }
      this_->rhsQS(Ns, t, x, xF, qdot, qdotF, tmp1, tmp2);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhsQS failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::rhsB(double t, const double* x, const double *rx, double* rxdot) {
    log("CvodesInterface::rhsB", "begin");

    // Pass inputs
    g_.setInputNZ(&t, RDAE_T);
    g_.setInputNZ(x, RDAE_X);
    g_.setInput(p(), RDAE_P);
    g_.setInput(rp(), RDAE_RP);
    g_.setInputNZ(rx, RDAE_RX);

    if (monitor_rhsB_) {
      userOut() << "t       = " << t << endl;
      userOut() << "x       = " << g_.input(RDAE_X) << endl;
      userOut() << "p       = " << g_.input(RDAE_P) << endl;
      userOut() << "rx      = " << g_.input(RDAE_RX) << endl;
      userOut() << "rp      = " << g_.input(RDAE_RP) << endl;
    }

    // Evaluate
    g_.evaluate();

    // Save to output
    g_.getOutputNZ(rxdot, RDAE_ODE);

    if (monitor_rhsB_) {
      userOut() << "xdotB = " << g_.output(RDAE_ODE) << endl;
    }

    // Negate (note definition of g)
    for (int i=0; i<nrx_; ++i)
      rxdot[i] *= -1;

    log("CvodesInterface::rhsB", "end");
  }

  void CvodesInterface::rhsBS(double t, N_Vector x, N_Vector *xF, N_Vector rx, N_Vector rxdot) {

    // Commented out since a new implementation currently cannot be tested
    casadi_error("Commented out, #884, #794.");
  }

  int CvodesInterface::rhsB_wrapper(double t, N_Vector x, N_Vector rx, N_Vector rxdot,
                                   void *user_data) {
    try {
      casadi_assert(user_data);
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      this_->rhsB(t, NV_DATA_S(x), NV_DATA_S(rx), NV_DATA_S(rxdot));
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhsB failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::rhsBS_wrapper(double t, N_Vector x, N_Vector *xF, N_Vector xB,
                                    N_Vector xdotB, void *user_data) {
    try {
      casadi_assert(user_data);
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      this_->rhsBS(t, x, xF, xB, xdotB);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhsBS failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::rhsQB_wrapper(double t, N_Vector x, N_Vector rx,
                                    N_Vector rqdot, void *user_data) {
    try {
      casadi_assert(user_data);
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      this_->rhsQB(t, NV_DATA_S(x), NV_DATA_S(rx), NV_DATA_S(rqdot));
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhsQB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::rhsQB(double t, const double* x, const double* rx, double* rqdot) {
    if (monitor_rhsQB_) {
      userOut() << "CvodesInterface::rhsQB: begin" << endl;
    }

    // Pass inputs
    g_.setInputNZ(&t, RDAE_T);
    g_.setInputNZ(x, RDAE_X);
    g_.setInput(p(), RDAE_P);
    g_.setInput(rp(), RDAE_RP);
    g_.setInputNZ(rx, RDAE_RX);

    if (monitor_rhsB_) {
      userOut() << "t       = " << t << endl;
      userOut() << "x       = " << g_.input(RDAE_X) << endl;
      userOut() << "p       = " << g_.input(RDAE_P) << endl;
      userOut() << "rx      = " << g_.input(RDAE_RX) << endl;
      userOut() << "rp      = " << g_.input(RDAE_RP) << endl;
    }

    // Evaluate
    g_.evaluate();

    // Save to output
    g_.getOutputNZ(rqdot, RDAE_QUAD);

    if (monitor_rhsB_) {
      userOut() << "qdotB = " << g_.output(RDAE_QUAD) << endl;
    }

    // Negate (note definition of g)
    for (int i=0; i<nrq_; ++i)
      rqdot[i] *= -1;

    if (monitor_rhsQB_) {
      userOut() << "CvodesInterface::rhsQB: end" << endl;
    }
  }

  int CvodesInterface::jtimes_wrapper(N_Vector v, N_Vector Jv, double t, N_Vector x,
                                     N_Vector xdot, void *user_data, N_Vector tmp) {
    try {
      casadi_assert(user_data);
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      this_->jtimes(v, Jv, t, x, xdot, tmp);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "jtimes failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::jtimesB_wrapper(N_Vector vB, N_Vector JvB, double t, N_Vector x,
                                      N_Vector xB, N_Vector xdotB, void *user_data ,
                                      N_Vector tmpB) {
    try {
      casadi_assert(user_data);
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      this_->jtimesB(vB, JvB, t, x, xB, xdotB, tmpB);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "jtimes failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::jtimes(N_Vector v, N_Vector Jv, double t, N_Vector x,
                              N_Vector xdot, N_Vector tmp) {
    log("CvodesInterface::jtimes", "begin");
    // Get time
    time1 = clock();

    // Pass input
    f_fwd_.setInputNZ(&t,                 DAE_T);
    f_fwd_.setInputNZ(NV_DATA_S(x),       DAE_X);
    f_fwd_.setInput(p(), DAE_P);

    // Pass input seeds
    f_fwd_.setInput(0.0,          DAE_NUM_IN + DAE_T);
    f_fwd_.setInputNZ(NV_DATA_S(v), DAE_NUM_IN + DAE_X);
    f_fwd_.setInput(0.0,          DAE_NUM_IN + DAE_P);

    // Evaluate
    f_fwd_.evaluate();

    // Get the output seeds
    f_fwd_.getOutputNZ(NV_DATA_S(Jv), DAE_NUM_OUT + DAE_ODE);

    // Log time duration
    time2 = clock();
    t_jac += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;

    log("CvodesInterface::jtimes", "end");
  }

  void CvodesInterface::jtimesB(N_Vector vB, N_Vector JvB, double t, N_Vector x, N_Vector xB,
                               N_Vector xdotB, N_Vector tmpB) {
    log("CvodesInterface::jtimesB", "begin");
    // Get time
    time1 = clock();

    // Pass input
    g_fwd_.setInputNZ(&t,                  RDAE_T);
    g_fwd_.setInputNZ(NV_DATA_S(x),        RDAE_X);
    g_fwd_.setInput(p(), RDAE_P);
    g_fwd_.setInputNZ(NV_DATA_S(xB),       RDAE_RX);
    g_fwd_.setInput(rp(), RDAE_RP);

    // Pass input seeds
    g_fwd_.setInput(0.0,           RDAE_NUM_IN + RDAE_T);
    g_fwd_.setInput(0.0,           RDAE_NUM_IN + RDAE_X);
    g_fwd_.setInput(0.0,           RDAE_NUM_IN + RDAE_P);
    g_fwd_.setInputNZ(NV_DATA_S(vB), RDAE_NUM_IN + RDAE_RX);
    g_fwd_.setInput(0.0,           RDAE_NUM_IN + RDAE_RP);

    // Evaluate
    g_fwd_.evaluate();

    // Get the output seeds
    g_fwd_.getOutputNZ(NV_DATA_S(JvB), RDAE_NUM_OUT + RDAE_ODE);

    // Log time duration
    time2 = clock();
    t_jac += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
    log("CvodesInterface::jtimesB", "end");
  }

  int CvodesInterface::djac_wrapper(long N, double t, N_Vector x, N_Vector xdot, DlsMat Jac,
                                   void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    try {
      casadi_assert(user_data);
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      this_->djac(N, t, x, xdot, Jac, tmp1, tmp2, tmp3);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "djac failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::djacB_wrapper(long NeqB, double t, N_Vector x, N_Vector xB, N_Vector xdotB,
                                    DlsMat JacB, void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                                    N_Vector tmp3B) {
    try {
      casadi_assert(user_data);
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      this_->djacB(NeqB, t, x, xB, xdotB, JacB, tmp1B, tmp2B, tmp3B);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "djacB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::djac(long N, double t, N_Vector x, N_Vector xdot, DlsMat Jac,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    log("CvodesInterface::djac", "begin");

    // Get time
    time1 = clock();

    // Pass inputs to the jacobian function
    jac_.setInputNZ(&t, DAE_T);
    jac_.setInputNZ(NV_DATA_S(x), DAE_X);
    jac_.setInput(f_.input(DAE_P), DAE_P);
    jac_.setInput(1.0, DAE_NUM_IN);
    jac_.setInput(0.0, DAE_NUM_IN+1);


    if (monitored("djac")) {
      userOut() << "DAE_T    = " << t << endl;
      userOut() << "DAE_X    = " << jac_.input(DAE_X) << endl;
      userOut() << "RDAE_P    = " << jac_.input(DAE_P) << endl;
      userOut() << "c_x = " << jac_.input(DAE_NUM_IN) << endl;
      userOut() << "c_xdot = " << jac_.input(DAE_NUM_IN+1) << endl;
    }

    // Evaluate
    jac_.evaluate();

    if (monitored("djac")) {
      userOut() << "jac = " << jac_.output() << endl;
    }

    // Get sparsity and non-zero elements
    const int* colind = jac_.sparsity_out(0).colind();
    int ncol = jac_.size2_out(0);
    const int* row = jac_.sparsity_out(0).row();
    const vector<double>& val = jac_.output().data();

    // Loop over columns
    for (int cc=0; cc<ncol; ++cc) {
      // Loop over non-zero entries
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        // Get row
        int rr = row[el];

        // Set the element
        DENSE_ELEM(Jac, rr, cc) = val[el];
      }
    }

    // Log time duration
    time2 = clock();
    t_jac += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;

    log("CvodesInterface::djac", "end");
  }

  void CvodesInterface::djacB(long NeqB, double t, N_Vector x, N_Vector xB, N_Vector xdotB,
                             DlsMat JacB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
    log("CvodesInterface::djacB", "begin");
    // Get time
    time1 = clock();

    // Pass inputs to the jacobian function
    jacB_.setInputNZ(&t, RDAE_T);
    jacB_.setInputNZ(NV_DATA_S(x), RDAE_X);
    jacB_.setInput(p(), RDAE_P);
    jacB_.setInputNZ(NV_DATA_S(xB), RDAE_RX);
    jacB_.setInput(rp(), RDAE_RP);
    jacB_.setInput(-1.0, RDAE_NUM_IN); // validated
    jacB_.setInput(0.0, RDAE_NUM_IN+1);


    if (monitored("djacB")) {
      userOut() << "RDAE_T    = " << t << endl;
      userOut() << "RDAE_X    = " << jacB_.input(RDAE_X) << endl;
      userOut() << "RDAE_P    = " << jacB_.input(RDAE_P) << endl;
      userOut() << "RDAE_RX    = " << jacB_.input(RDAE_RX) << endl;
      userOut() << "RDAE_RP    = " << jacB_.input(RDAE_RP) << endl;
      userOut() << "c_x = " << jacB_.input(RDAE_NUM_IN) << endl;
      userOut() << "c_xdot = " << jacB_.input(RDAE_NUM_IN+1) << endl;
    }

    // Evaluate
    jacB_.evaluate();

    if (monitored("djacB")) {
      userOut() << "jacB = " << jacB_.output() << endl;
    }

    // Get sparsity and non-zero elements
    const int* colind = jacB_.sparsity_out(0).colind();
    int ncol = jacB_.size2_out(0);
    const int* row = jacB_.sparsity_out(0).row();
    const vector<double>& val = jacB_.output().data();

    // Loop over columns
    for (int cc=0; cc<ncol; ++cc) {
      // Loop over non-zero entries
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        // Get row
        int rr = row[el];

        // Set the element
        DENSE_ELEM(JacB, rr, cc) = val[el];
      }
    }

    // Log time duration
    time2 = clock();
    t_jac += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
    log("CvodesInterface::djacB", "end");
  }

  int CvodesInterface::bjac_wrapper(long N, long mupper, long mlower, double t, N_Vector x,
                                   N_Vector xdot, DlsMat Jac, void *user_data,
                                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    try {
      casadi_assert(user_data);
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      this_->bjac(N, mupper, mlower, t, x, xdot, Jac, tmp1, tmp2, tmp3);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "bjac failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::bjacB_wrapper(long NeqB, long mupperB, long mlowerB, double t, N_Vector x,
                                    N_Vector xB, N_Vector xdotB, DlsMat JacB, void *user_data,
                                    N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
    try {
      casadi_assert(user_data);
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      this_->bjacB(NeqB, mupperB, mlowerB, t, x, xB, xdotB, JacB, tmp1B, tmp2B, tmp3B);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "bjacB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::bjac(long N, long mupper, long mlower, double t, N_Vector x, N_Vector xdot,
                            DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    log("CvodesInterface::bjac", "begin");

    // Get time
    time1 = clock();

    // Pass inputs to the jacobian function
    jac_.setInputNZ(&t, DAE_T);
    jac_.setInputNZ(NV_DATA_S(x), DAE_X);
    jac_.setInput(f_.input(DAE_P), DAE_P);
    jac_.setInput(1.0, DAE_NUM_IN);
    jac_.setInput(0.0, DAE_NUM_IN+1);

    // Evaluate
    jac_.evaluate();

    // Get sparsity and non-zero elements
    const int* colind = jac_.sparsity_out(0).colind();
    int ncol = jac_.size2_out(0);
    const int* row = jac_.sparsity_out(0).row();
    const vector<double>& val = jac_.output().data();

    // Loop over cols
    for (int cc=0; cc<ncol; ++cc) {
      // Loop over non-zero entries
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        // Get row
        int rr = row[el];

        // Set the element
        if (cc-rr<=mupper && rr-cc<=mlower)
          BAND_ELEM(Jac, rr, cc) = val[el];
      }
    }

    // Log time duration
    time2 = clock();
    t_jac += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;

    log("CvodesInterface::bjac", "end");
  }

  void CvodesInterface::bjacB(long NeqB, long mupperB, long mlowerB, double t, N_Vector x,
                             N_Vector xB, N_Vector xdotB, DlsMat JacB, N_Vector tmp1B,
                             N_Vector tmp2B, N_Vector tmp3B) {
    log("CvodesInterface::bjacB", "begin");

    // Get time
    time1 = clock();

    // Pass inputs to the jacobian function
    jacB_.setInputNZ(&t, RDAE_T);
    jacB_.setInputNZ(NV_DATA_S(x), RDAE_X);
    jacB_.setInput(p(), RDAE_P);
    jacB_.setInputNZ(NV_DATA_S(xB), RDAE_RX);
    jacB_.setInput(rp(), RDAE_RP);
    jacB_.setInput(-1.0, RDAE_NUM_IN);
    jacB_.setInput(0.0, RDAE_NUM_IN+1);

    if (monitored("bjacB")) {
      userOut() << "RDAE_T    = " << t << endl;
      userOut() << "RDAE_X    = " << jacB_.input(RDAE_X) << endl;
      userOut() << "RDAE_P    = " << jacB_.input(RDAE_P) << endl;
      userOut() << "RDAE_RX    = " << jacB_.input(RDAE_RX) << endl;
      userOut() << "RDAE_RP    = " << jacB_.input(RDAE_RP) << endl;
      userOut() << "c_x = " << jacB_.input(RDAE_NUM_IN) << endl;
      userOut() << "c_xdot = " << jacB_.input(RDAE_NUM_IN+1) << endl;
    }

    // Evaluate
    jacB_.evaluate();

    if (monitored("bjacB")) {
      userOut() << "jacB = " << jacB_.output() << endl;
    }

    // Get sparsity and non-zero elements
    const int* colind = jacB_.sparsity_out(0).colind();
    int ncol = jacB_.size2_out(0);
    const int* row = jacB_.sparsity_out(0).row();
    const vector<double>& val = jacB_.output().data();

    // Loop over columns
    for (int cc=0; cc<ncol; ++cc) {
      // Loop over non-zero entries
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        // Get row
        int rr = row[el];

        // Set the element
        if (cc-rr<=mupperB && rr-cc<=mlowerB)
          BAND_ELEM(JacB, rr, cc) = val[el];
      }
    }

    // Log time duration
    time2 = clock();
    t_jac += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;

    log("CvodesInterface::bjacB", "end");
  }

  void CvodesInterface::setStopTime(double tf) {
    // Set the stop time of the integration -- don't integrate past this point
    int flag = CVodeSetStopTime(mem_, tf);
    if (flag != CV_SUCCESS) cvodes_error("CVodeSetStopTime", flag);
  }

  int CvodesInterface::psolve_wrapper(double t, N_Vector x, N_Vector xdot, N_Vector r,
                                     N_Vector z, double gamma, double delta, int lr,
                                     void *user_data, N_Vector tmp) {
    try {
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      casadi_assert(this_);
      this_->psolve(t, x, xdot, r, z, gamma, delta, lr, tmp);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "psolve failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::psolveB_wrapper(double t, N_Vector x, N_Vector xB, N_Vector xdotB,
                                      N_Vector rvecB, N_Vector zvecB, double gammaB,
                                      double deltaB, int lr, void *user_data, N_Vector tmpB) {
    try {
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      casadi_assert(this_);
      this_->psolveB(t, x, xB, xdotB, rvecB, zvecB, gammaB, deltaB, lr, tmpB);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "psolveB failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::psetup_wrapper(double t, N_Vector x, N_Vector xdot, booleantype jok,
                                     booleantype *jcurPtr, double gamma, void *user_data,
                                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    try {
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      casadi_assert(this_);
      this_->psetup(t, x, xdot, jok, jcurPtr, gamma, tmp1, tmp2, tmp3);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "psetup failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::psetupB_wrapper(double t, N_Vector x, N_Vector xB, N_Vector xdotB,
                                      booleantype jokB, booleantype *jcurPtrB, double gammaB,
                                      void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                                      N_Vector tmp3B) {
    try {
      CvodesInterface *this_ = static_cast<CvodesInterface*>(user_data);
      casadi_assert(this_);
      this_->psetupB(t, x, xB, xdotB, jokB, jcurPtrB, gammaB, tmp1B, tmp2B, tmp3B);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "psetupB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::psolve(double t, N_Vector x, N_Vector xdot, N_Vector r, N_Vector z,
                              double gamma, double delta, int lr, N_Vector tmp) {
    // Get time
    time1 = clock();

    // Copy input to output, if necessary
    if (r!=z) {
      N_VScale(1.0, r, z);
    }

    // Solve the (possibly factorized) system
    casadi_assert(linsol_.nnz_out(0) == NV_LENGTH_S(z));
    linsol_.solve(NV_DATA_S(z), 1, false);

    // Log time duration
    time2 = clock();
    t_lsolve += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
  }

  void CvodesInterface::psolveB(double t, N_Vector x, N_Vector xB, N_Vector xdotB, N_Vector rvecB,
                               N_Vector zvecB, double gammaB, double deltaB,
                               int lr, N_Vector tmpB) {
    // Get time
    time1 = clock();

    // Copy input to output, if necessary
    if (rvecB!=zvecB) {
      N_VScale(1.0, rvecB, zvecB);
    }

    // Solve the (possibly factorized) system
    casadi_assert(linsolB_.nnz_out(0) == NV_LENGTH_S(zvecB));
    linsolB_.solve(NV_DATA_S(zvecB), 1, false);

    // Log time duration
    time2 = clock();
    t_lsolve += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
  }

  void CvodesInterface::psetup(double t, N_Vector x, N_Vector xdot, booleantype jok,
                              booleantype *jcurPtr, double gamma, N_Vector tmp1,
                              N_Vector tmp2, N_Vector tmp3) {
    log("CvodesInterface::psetup", "begin");
    // Get time
    time1 = clock();

    // Pass input to the jacobian function
    jac_.setInputNZ(&t, DAE_T);
    jac_.setInputNZ(NV_DATA_S(x), DAE_X);
    jac_.setInput(p(), DAE_P);
    jac_.setInput(-gamma, DAE_NUM_IN);
    jac_.setInput(1.0, DAE_NUM_IN+1);

    // Evaluate jacobian
    jac_.evaluate();

    // Log time duration
    time2 = clock();
    t_lsetup_jac += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;

    // Pass non-zero elements, scaled by -gamma, to the linear solver
    linsol_.setInput(jac_.output(), 0);

    // Prepare the solution of the linear system (e.g. factorize)
    // -- only if the linear solver inherits from LinearSolver
    linsol_.prepare();

    // Log time duration
    time1 = clock();
    t_lsetup_fac += static_cast<double>(time1-time2)/CLOCKS_PER_SEC;

    log("CvodesInterface::psetup", "end");
  }

  void CvodesInterface::psetupB(double t, N_Vector x, N_Vector xB, N_Vector xdotB,
                               booleantype jokB, booleantype *jcurPtrB, double gammaB,
                               N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
    log("CvodesInterface::psetupB", "begin");
    // Get time
    time1 = clock();

    // Pass inputs to the jacobian function
    jacB_.setInputNZ(&t, RDAE_T);
    jacB_.setInputNZ(NV_DATA_S(x), RDAE_X);
    jacB_.setInput(p(), RDAE_P);
    jacB_.setInputNZ(NV_DATA_S(xB), RDAE_RX);
    jacB_.setInput(rp(), RDAE_RP);
    jacB_.setInput(gammaB, RDAE_NUM_IN); // validated
    jacB_.setInput(1.0, RDAE_NUM_IN+1); // validated

    if (monitored("psetupB")) {
      userOut() << "RDAE_T    = " << t << endl;
      userOut() << "RDAE_X    = " << jacB_.input(RDAE_X) << endl;
      userOut() << "RDAE_P    = " << jacB_.input(RDAE_P) << endl;
      userOut() << "RDAE_RX    = " << jacB_.input(RDAE_RX) << endl;
      userOut() << "RDAE_RP    = " << jacB_.input(RDAE_RP) << endl;
      userOut() << "gamma = " << gammaB << endl;
    }

    // Evaluate jacobian
    jacB_.evaluate();

    if (monitored("psetupB")) {
      userOut() << "psetupB = " << jacB_.output() << endl;
    }

    // Log time duration
    time2 = clock();
    t_lsetup_jac += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;

    // Pass non-zero elements, scaled by -gamma, to the linear solver
    linsolB_.setInput(jacB_.output(), 0);

    // Prepare the solution of the linear system (e.g. factorize)
    // -- only if the linear solver inherits from LinearSolver
    linsolB_.prepare();

    // Log time duration
    time1 = clock();
    t_lsetup_fac += static_cast<double>(time1-time2)/CLOCKS_PER_SEC;
    log("CvodesInterface::psetupB", "end");
  }

  void CvodesInterface::lsetup(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
                              booleantype *jcurPtr,
                              N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3) {
    // Current time
    double t = cv_mem->cv_tn;

    // Scaling factor before J
    double gamma = cv_mem->cv_gamma;

    // Call the preconditioner setup function (which sets up the linear solver)
    psetup(t, x, xdot, FALSE, jcurPtr, gamma, vtemp1, vtemp2, vtemp3);
  }

  void CvodesInterface::lsetupB(double t, double gamma, int convfail,
                               N_Vector x, N_Vector xB, N_Vector xdotB, booleantype *jcurPtr,
                               N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3) {
    // Call the preconditioner setup function (which sets up the linear solver)
    psetupB(t, x, xB, xdotB, FALSE, jcurPtr, gamma, vtemp1, vtemp2, vtemp3);
  }

  int CvodesInterface::lsetup_wrapper(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
                                     booleantype *jcurPtr,
                                     N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3) {
    try {
      CvodesInterface *this_ = static_cast<CvodesInterface*>(cv_mem->cv_lmem);
      casadi_assert(this_);
      this_->lsetup(cv_mem, convfail, x, xdot, jcurPtr, vtemp1, vtemp2, vtemp3);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "lsetup failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::lsetupB_wrapper(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
                                      booleantype *jcurPtr,
                                      N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3) {
    try {
      CvodesInterface *this_ = static_cast<CvodesInterface*>(cv_mem->cv_lmem);
      casadi_assert(this_);
      CVadjMem ca_mem;
      //CVodeBMem cvB_mem;

      int flag;

      // Current time
      double t = cv_mem->cv_tn; // TODO(Joel): is this correct?
      double gamma = cv_mem->cv_gamma;

      cv_mem = static_cast<CVodeMem>(cv_mem->cv_user_data);

      ca_mem = cv_mem->cv_adj_mem;
      //cvB_mem = ca_mem->ca_bckpbCrt;

      // Get FORWARD solution from interpolation.
      flag = ca_mem->ca_IMget(cv_mem, t, ca_mem->ca_ytmp, NULL);
      if (flag != CV_SUCCESS) casadi_error("Could not interpolate forward states");

      this_->lsetupB(t, gamma, convfail, ca_mem->ca_ytmp, x, xdot, jcurPtr, vtemp1, vtemp2, vtemp3);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "lsetupB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::lsolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                              N_Vector x, N_Vector xdot) {
    log("CvodesInterface::lsolve", "begin");

    // Current time
    double t = cv_mem->cv_tn;

    // Scaling factor before J
    double gamma = cv_mem->cv_gamma;

    // Accuracy
    double delta = 0.0;

    // Left/right preconditioner
    int lr = 1;

    // Call the preconditioner solve function (which solves the linear system)
    psolve(t, x, xdot, b, b, gamma, delta, lr, 0);

    log("CvodesInterface::lsolve", "end");
  }

  void CvodesInterface::lsolveB(double t, double gamma, N_Vector b, N_Vector weight,
                               N_Vector x, N_Vector xB, N_Vector xdotB) {
    log("CvodesInterface::lsolveB", "begin");
    // Accuracy
    double delta = 0.0;

    // Left/right preconditioner
    int lr = 1;

    // Call the preconditioner solve function (which solves the linear system)
    psolveB(t, x, xB, xdotB, b, b, gamma, delta, lr, 0);

    log("CvodesInterface::lsolveB", "end");
  }

  int CvodesInterface::lsolve_wrapper(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                                     N_Vector x, N_Vector xdot) {
    try {
      CvodesInterface *this_ = static_cast<CvodesInterface*>(cv_mem->cv_lmem);
      casadi_assert(this_);
      this_->lsolve(cv_mem, b, weight, x, xdot);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "lsolve failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::lsolveB_wrapper(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                                      N_Vector x, N_Vector xdot) {
    try {
      CvodesInterface *this_ = static_cast<CvodesInterface*>(cv_mem->cv_lmem);
      casadi_assert(this_);
      CVadjMem ca_mem;
      //CVodeBMem cvB_mem;

      int flag;

      // Current time
      double t = cv_mem->cv_tn; // TODO(Joel): is this correct?
      double gamma = cv_mem->cv_gamma;

      cv_mem = static_cast<CVodeMem>(cv_mem->cv_user_data);

      ca_mem = cv_mem->cv_adj_mem;
      //cvB_mem = ca_mem->ca_bckpbCrt;

      // Get FORWARD solution from interpolation.
      flag = ca_mem->ca_IMget(cv_mem, t, ca_mem->ca_ytmp, NULL);
      if (flag != CV_SUCCESS) casadi_error("Could not interpolate forward states");

      this_->lsolveB(t, gamma, b, weight, ca_mem->ca_ytmp, x, xdot);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "lsolveB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::initDenseLinearSolver() {
    int flag = CVDense(mem_, nx_);
    if (flag!=CV_SUCCESS) cvodes_error("CVDense", flag);
    if (exact_jacobian_) {
      flag = CVDlsSetDenseJacFn(mem_, djac_wrapper);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsSetDenseJacFn", flag);
    }
  }

  void CvodesInterface::initBandedLinearSolver() {
    pair<int, int> bw = getBandwidth();
    int flag = CVBand(mem_, nx_, bw.first, bw.second);
    if (flag!=CV_SUCCESS) cvodes_error("CVBand", flag);
    if (exact_jacobian_) {
      flag = CVDlsSetBandJacFn(mem_, bjac_wrapper);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsSetBandJacFn", flag);
    }
  }

  void CvodesInterface::initIterativeLinearSolver() {
    // Attach the sparse solver
    int flag;
    switch (itsol_f_) {
    case SD_GMRES:
      flag = CVSpgmr(mem_, pretype_f_, max_krylov_);
      if (flag!=CV_SUCCESS) cvodes_error("CVBand", flag);
      break;
    case SD_BCGSTAB:
      flag = CVSpbcg(mem_, pretype_f_, max_krylov_);
      if (flag!=CV_SUCCESS) cvodes_error("CVSpbcg", flag);
      break;
    case SD_TFQMR:
      flag = CVSptfqmr(mem_, pretype_f_, max_krylov_);
      if (flag!=CV_SUCCESS) cvodes_error("CVSptfqmr", flag);
      break;
    }

    // Attach functions for jacobian information
    if (exact_jacobian_) {
      // Form the Jacobian-times-vector function
      f_fwd_ = f_.derivative(1, 0);

      flag = CVSpilsSetJacTimesVecFn(mem_, jtimes_wrapper);
      if (flag!=CV_SUCCESS) cvodes_error("CVSpilsSetJacTimesVecFn", flag);
    }

    // Add a preconditioner
    if (use_preconditioner_) {
      // Make sure that a Jacobian has been provided
      if (jac_.isNull())
          throw CasadiException("CvodesInterface::init(): No Jacobian has been provided.");

      // Make sure that a linear solver has been provided
      if (linsol_.isNull())
          throw CasadiException("CvodesInterface::init(): "
                                "No user defined linear solver has been provided.");

      // Pass to IDA
      flag = CVSpilsSetPreconditioner(mem_, psetup_wrapper, psolve_wrapper);
      if (flag != CV_SUCCESS) cvodes_error("CVSpilsSetPreconditioner", flag);
    }
  }

  void CvodesInterface::initUserDefinedLinearSolver() {
    // Make sure that a Jacobian has been provided
    if (jac_.isNull())
        throw CasadiException("CvodesInterface::initUserDefinedLinearSolver(): "
                              "No Jacobian has been provided.");

    // Make sure that a linear solver has been provided
    if (linsol_.isNull())
        throw CasadiException("CvodesInterface::initUserDefinedLinearSolver(): "
                              "No user defined linear solver has been provided.");

    //  Set fields in the IDA memory
    CVodeMem cv_mem = static_cast<CVodeMem>(mem_);
    cv_mem->cv_lmem   = this;
    cv_mem->cv_lsetup = lsetup_wrapper;
    cv_mem->cv_lsolve = lsolve_wrapper;
    cv_mem->cv_setupNonNull = TRUE;
  }

  void CvodesInterface::initDenseLinearSolverB() {
    int flag = CVDenseB(mem_, whichB_, nrx_);
    if (flag!=CV_SUCCESS) cvodes_error("CVDenseB", flag);
    if (exact_jacobianB_) {
      flag = CVDlsSetDenseJacFnB(mem_, whichB_, djacB_wrapper);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsSetDenseJacFnB", flag);
    }
  }

  void CvodesInterface::initBandedLinearSolverB() {
    pair<int, int> bw = getBandwidthB();
    int flag = CVBandB(mem_, whichB_, nrx_, bw.first, bw.second);
    if (flag!=CV_SUCCESS) cvodes_error("CVBandB", flag);

    if (exact_jacobianB_) {
      flag = CVDlsSetBandJacFnB(mem_, whichB_, bjacB_wrapper);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsSetBandJacFnB", flag);
    }
  }

  void CvodesInterface::initIterativeLinearSolverB() {
    int flag;
    switch (itsol_g_) {
    case SD_GMRES:
      flag = CVSpgmrB(mem_, whichB_, pretype_g_, max_krylovB_);
      if (flag!=CV_SUCCESS) cvodes_error("CVSpgmrB", flag);
      break;
    case SD_BCGSTAB:
      flag = CVSpbcgB(mem_, whichB_, pretype_g_, max_krylovB_);
      if (flag!=CV_SUCCESS) cvodes_error("CVSpbcgB", flag);
      break;
    case SD_TFQMR:
      flag = CVSptfqmrB(mem_, whichB_, pretype_g_, max_krylovB_);
      if (flag!=CV_SUCCESS) cvodes_error("CVSptfqmrB", flag);
      break;
    }

    // Attach functions for jacobian information
    if (exact_jacobianB_) {
      // Form the Jacobian-times-vector function
      g_fwd_ = g_.derivative(1, 0);

      flag = CVSpilsSetJacTimesVecFnB(mem_, whichB_, jtimesB_wrapper);
      if (flag!=CV_SUCCESS) cvodes_error("CVSpilsSetJacTimesVecFnB", flag);
    }

    // Add a preconditioner
    if (use_preconditionerB_) {
      // Make sure that a Jacobian has been provided
      if (jacB_.isNull())
        casadi_error("CvodesInterface::init(): No backwards Jacobian has been provided.");

      // Make sure that a linear solver has been provided
      if (linsolB_.isNull())
        casadi_error("CvodesInterface::init(): "
                     "No user defined backwards  linear solver has been provided.");

      // Pass to IDA
      flag = CVSpilsSetPreconditionerB(mem_, whichB_, psetupB_wrapper, psolveB_wrapper);
      if (flag != CV_SUCCESS) cvodes_error("CVSpilsSetPreconditionerB", flag);
    }

  }

  void CvodesInterface::initUserDefinedLinearSolverB() {
    // Make sure that a Jacobian has been provided
    if (jacB_.isNull())
        throw CasadiException("CvodesInterface::initUserDefinedLinearSolverB(): "
                              "No backwards Jacobian has been provided.");

    // Make sure that a linear solver has been provided
    if (linsolB_.isNull())
        throw CasadiException("CvodesInterface::initUserDefinedLinearSolverB(): "
                              "No user defined backward linear solver has been provided.");

    CVodeMem cv_mem = static_cast<CVodeMem>(mem_);
    CVadjMem ca_mem = cv_mem->cv_adj_mem;
    CVodeBMem cvB_mem = ca_mem->cvB_mem;
    cvB_mem->cv_lmem   = this;

    cvB_mem->cv_mem->cv_lmem = this;
    cvB_mem->cv_mem->cv_lsetup = lsetupB_wrapper;
    cvB_mem->cv_mem->cv_lsolve = lsolveB_wrapper;
    cvB_mem->cv_mem->cv_setupNonNull = TRUE;
  }

  template<typename MatType>
  Function CvodesInterface::getJacGen() {
    // Get the Jacobian in the Newton iteration
    MatType c_x = MatType::sym("c_x");
    MatType c_xdot = MatType::sym("c_xdot");
    MatType jac = c_x*MatType::jac(f_, DAE_X, DAE_ODE) + c_xdot*MatType::eye(nx_);

    // Jacobian function
    std::vector<MatType> jac_in = MatType::get_input(f_);
    jac_in.push_back(c_x);
    jac_in.push_back(c_xdot);

    // Return generated function
    return Function("jac", jac_in, {jac});
  }

  template<typename MatType>
  Function CvodesInterface::getJacGenB() {
    // Get the Jacobian in the Newton iteration
    MatType c_x = MatType::sym("c_x");
    MatType c_xdot = MatType::sym("c_xdot");
    MatType jac = c_x*MatType::jac(g_, RDAE_RX, RDAE_ODE) + c_xdot*MatType::eye(nrx_);

    // Jacobian function
    std::vector<MatType> jac_in = MatType::get_input(g_);
    jac_in.push_back(c_x);
    jac_in.push_back(c_xdot);

    // return generated function
    return Function("jacB", jac_in, {jac});
  }

  Function CvodesInterface::getJacB() {
    if (g_.is_a("sxfunction")) {
      return getJacGenB<SX>();
    } else if (g_.is_a("mxfunction")) {
      return getJacGenB<MX>();
    } else {
      throw CasadiException("CvodesInterface::getJacB(): Not an SXFunction or MXFunction");
    }
  }


  Function CvodesInterface::getJac() {
    if (f_.is_a("sxfunction")) {
      return getJacGen<SX>();
    } else if (f_.is_a("mxfunction")) {
      return getJacGen<MX>();
    } else {
      throw CasadiException("CvodesInterface::getJac(): Not an SXFunction or MXFunction");
    }
  }

  CvodesInterface::Memory::Memory(CvodesInterface& s) : self(s) {
  }

  CvodesInterface::Memory::~Memory() {
  }



} // namespace casadi
