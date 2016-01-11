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

    disable_internal_warnings_ = false;
  }

  void CvodesInterface::freeCVodes() {
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

    // Algebraic variables not supported
    casadi_assert_message(nz_==0 && nrz_==0,
                          "CVODES does not support algebraic variables");

    // Read options
    monitor_rhsB_  = monitored("resB");
    monitor_rhs_   = monitored("res");
    monitor_rhsQB_ = monitored("resQB");

    if (option("linear_multistep_method")=="adams")  lmm_ = CV_ADAMS;
    else if (option("linear_multistep_method")=="bdf") lmm_ = CV_BDF;
    else
      throw CasadiException("Unknown linear multistep method");

    if (option("nonlinear_solver_iteration")=="newton") iter_ = CV_NEWTON;
    else if (option("nonlinear_solver_iteration")=="functional") iter_ = CV_FUNCTIONAL;
    else
      throw CasadiException("Unknown nonlinear solver iteration");

    // Disable internal warning messages?
    disable_internal_warnings_ = option("disable_internal_warnings");

    // Attach functions for jacobian information
    if (exact_jacobian_) {
      switch (linsol_f_) {
      case SD_ITERATIVE:
        f_fwd_ = f_.derivative(1, 0);
        alloc(f_fwd_);
        break;
      default: break;
      }
    }

    if (exact_jacobianB_) {
      switch (linsol_g_) {
      case SD_ITERATIVE:
        g_fwd_ = g_.derivative(1, 0);
        alloc(g_fwd_);
        break;
      default: break;
      }
    }
  }

  void CvodesInterface::initAdj(CvodesMemory& m) const {

    // Create backward problem (use the same lmm and iter)
    int flag = CVodeCreateB(m.mem, lmm_, iter_, &m.whichB);
    if (flag != CV_SUCCESS) cvodes_error("CVodeCreateB", flag);

    // Initialize the backward problem
    double tB0 = grid_.back();
    flag = CVodeInitB(m.mem, m.whichB, rhsB_wrapper, tB0, m.rxz);
    if (flag != CV_SUCCESS) cvodes_error("CVodeInitB", flag);

    // Set tolerances
    flag = CVodeSStolerancesB(m.mem, m.whichB, reltolB_, abstolB_);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeSStolerancesB", flag);

    // User data
    flag = CVodeSetUserDataB(m.mem, m.whichB, &m);
    if (flag != CV_SUCCESS) cvodes_error("CVodeSetUserDataB", flag);

    // attach linear solver to backward problem
    switch (linsol_g_) {
    case SD_DENSE:
      initDenseLinsolB(m);
      break;
    case SD_BANDED:
      initBandedLinsolB(m);
      break;
    case SD_ITERATIVE:
      initIterativeLinsolB(m);
      break;
    case SD_USER_DEFINED:
      initUserDefinedLinsolB(m);
      break;
    }

    // Quadratures for the backward problem
    flag = CVodeQuadInitB(m.mem, m.whichB, rhsQB_wrapper, m.rq);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeQuadInitB", flag);

    if (option("quad_err_con").to_int()) {
      flag = CVodeSetQuadErrConB(m.mem, m.whichB, true);
      if (flag != CV_SUCCESS) cvodes_error("CVodeSetQuadErrConB", flag);

      flag = CVodeQuadSStolerancesB(m.mem, m.whichB, reltolB_, abstolB_);
      if (flag != CV_SUCCESS) cvodes_error("CVodeQuadSStolerancesB", flag);
    }

    // Mark initialized
    m.isInitAdj = true;
  }

  void CvodesInterface::init_memory(Memory& mem) const {
    SundialsInterface::init_memory(mem);
    CvodesMemory& m = dynamic_cast<CvodesMemory&>(mem);

    // Create CVodes memory block
    m.mem = CVodeCreate(lmm_, iter_);
    casadi_assert_message(m.mem!=0, "CVodeCreate: Creation failed");

    // Set error handler function
    int flag = CVodeSetErrHandlerFn(m.mem, ehfun_wrapper, &m);
    if (flag != CV_SUCCESS) cvodes_error("CVodeSetErrHandlerFn", flag);

    // Initialize CVodes
    double t0 = 0;
    flag = CVodeInit(m.mem, rhs_wrapper, t0, m.xz);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeInit", flag);

    // Set tolerances
    flag = CVodeSStolerances(m.mem, reltol_, abstol_);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeInit", flag);

    // Maximum number of steps
    CVodeSetMaxNumSteps(m.mem, option("max_num_steps").to_int());
    if (flag != CV_SUCCESS) cvodes_error("CVodeSetMaxNumSteps", flag);

    // attach a linear solver
    switch (linsol_f_) {
    case SD_DENSE:
      initDenseLinsol(m);
      break;
    case SD_BANDED:
      initBandedLinsol(m);
      break;
    case SD_ITERATIVE:
      initIterativeLinsol(m);
      break;
    case SD_USER_DEFINED:
      initUserDefinedLinsol(m);
      break;
    }

    // Set user data
    flag = CVodeSetUserData(m.mem, &m);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeSetUserData", flag);

    // Quadrature equations
    if (nq_>0) {
      // Initialize quadratures in CVodes
      flag = CVodeQuadInit(m.mem, rhsQ_wrapper, m.q);
      if (flag != CV_SUCCESS) cvodes_error("CVodeQuadInit", flag);

      // Should the quadrature errors be used for step size control?
      if (option("quad_err_con").to_int()) {
        flag = CVodeSetQuadErrCon(m.mem, true);
        if (flag != CV_SUCCESS) cvodes_error("CVodeSetQuadErrCon", flag);

        // Quadrature error tolerances
        // TODO(Joel): vector absolute tolerances
        flag = CVodeQuadSStolerances(m.mem, reltol_, abstol_);
        if (flag != CV_SUCCESS) cvodes_error("CVodeQuadSStolerances", flag);
      }
    }

    // Adjoint sensitivity problem
    if (!g_.isNull()) {
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
      flag = CVodeAdjInit(m.mem, Nd, interpType);
      if (flag != CV_SUCCESS) cvodes_error("CVodeAdjInit", flag);
      m.isInitAdj = false;
    }
  }

  void CvodesInterface::rhs(CvodesMemory& m, double t, N_Vector x, N_Vector xdot) const {
    log("CvodesInterface::rhs", "begin");

    // Get time
    m.time1 = clock();

    // Debug output
    if (monitor_rhs_) {
      printvar("t", t);
      printvar("x", x);
    }

    // Evaluate f_
    m.arg[DAE_T] = &t;
    m.arg[DAE_X] = NV_DATA_S(x);
    m.arg[DAE_Z] = 0;
    m.arg[DAE_P] = get_ptr(m.p);
    m.res[DAE_ODE] = NV_DATA_S(xdot);
    m.res[DAE_ALG] = 0;
    m.res[DAE_QUAD] = 0;
    f_(m.arg, m.res, m.iw, m.w, 0);

    // Debug output
    if (monitor_rhs_) {
      printvar("xdot", xdot);
    }

    // Log time
    m.time2 = clock();
    m.t_res += static_cast<double>(m.time2-m.time1)/CLOCKS_PER_SEC;

    log("CvodesInterface::rhs", "end");

  }

  int CvodesInterface::rhs_wrapper(double t, N_Vector x, N_Vector xdot, void *user_data) {
    try {
      casadi_assert(user_data);
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      m->self.rhs(*m, t, x, xdot);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhs failed: " << e.what() << endl;
      return 1;
    }
  }

  void CvodesInterface::reset(IntegratorMemory& mem, double t, const double* x,
                              const double* z, const double* _p) const {
    casadi_msg("CvodesInterface::reset begin");
    CvodesMemory& m = dynamic_cast<CvodesMemory&>(mem);

    // Reset the base classes
    SundialsInterface::reset(mem, t, x, z, _p);

    // Reset timers
    m.t_res = m.t_fres = m.t_jac = m.t_lsolve = m.t_lsetup_jac = m.t_lsetup_fac = 0;

    // Re-initialize
    int flag = CVodeReInit(m.mem, t, m.xz);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeReInit", flag);

    // Re-initialize quadratures
    if (nq_>0) {
      N_VConst(0.0, m.q);
      flag = CVodeQuadReInit(m.mem, m.q);
      if (flag != CV_SUCCESS) cvodes_error("CVodeQuadReInit", flag);
    }

    // Turn off sensitivities
    flag = CVodeSensToggleOff(m.mem);
    if (flag != CV_SUCCESS) cvodes_error("CVodeSensToggleOff", flag);

    // Re-initialize backward integration
    if (nrx_>0) {
      flag = CVodeAdjReInit(m.mem);
      if (flag != CV_SUCCESS) cvodes_error("CVodeAdjReInit", flag);
    }

    // Set the stop time of the integration -- don't integrate past this point
    if (stop_at_end_) setStopTime(m, grid_.back());
    casadi_msg("CvodesInterface::reset end");
  }

  void CvodesInterface::advance(IntegratorMemory& mem, double t, double* x,
                                double* z, double* q) const {
    CvodesMemory& m = dynamic_cast<CvodesMemory&>(mem);

    casadi_assert_message(t>=grid_.front(),
                          "CvodesInterface::integrate(" << t << "): "
                          "Cannot integrate to a time earlier than t0 ("
                          << grid_.front() << ")");
    casadi_assert_message(t<=grid_.back() || !stop_at_end_, "CvodesInterface::integrate("
                          << t << "):"
                          " Cannot integrate past a time later than tf (" << grid_.back() << ") "
                          "unless stop_at_end is set to False.");

    // Integrate, unless already at desired time
    const double ttol = 1e-9;
    if (fabs(m.t-t)>=ttol) {
      // Integrate forward ...
      if (nrx_>0) {
        // ... with taping
        int flag = CVodeF(m.mem, t, m.xz, &m.t, CV_NORMAL, &m.ncheck);
        if (flag!=CV_SUCCESS && flag!=CV_TSTOP_RETURN) cvodes_error("CVodeF", flag);
      } else {
        // ... without taping
        int flag = CVode(m.mem, t, m.xz, &m.t, CV_NORMAL);
        if (flag!=CV_SUCCESS && flag!=CV_TSTOP_RETURN) cvodes_error("CVode", flag);
      }

      // Get quadratures
      if (nq_>0) {
        double tret;
        int flag = CVodeGetQuad(m.mem, &tret, m.q);
        if (flag!=CV_SUCCESS) cvodes_error("CVodeGetQuad", flag);
      }
    }

    // Set function outputs
    casadi_copy(NV_DATA_S(m.xz), nx_, x);
    casadi_copy(NV_DATA_S(m.q), nq_, q);

    // Print statistics
    if (option("print_stats")) printStats(m, userOut());

    if (gather_stats_) {
      int flag = CVodeGetIntegratorStats(m.mem, &m.nsteps, &m.nfevals, &m.nlinsetups,
                                         &m.netfails, &m.qlast, &m.qcur, &m.hinused,
                                         &m.hlast, &m.hcur, &m.tcur);
      if (flag!=CV_SUCCESS) cvodes_error("CVodeGetIntegratorStats", flag);
    }

    casadi_msg("CvodesInterface::integrate(" << t << ") end");
  }

  void CvodesInterface::resetB(IntegratorMemory& mem, double t, const double* rx,
                               const double* rz, const double* rp) const {
    CvodesMemory& m = dynamic_cast<CvodesMemory&>(mem);

    // Reset the base classes
    SundialsInterface::resetB(mem, t, rx, rz, rp);

    int flag;
    if (m.isInitAdj) {
      flag = CVodeReInitB(m.mem, m.whichB, grid_.back(), m.rxz);
      if (flag != CV_SUCCESS) cvodes_error("CVodeReInitB", flag);

      flag = CVodeQuadReInitB(m.mem, m.whichB, m.rq);
      if (flag!=CV_SUCCESS) cvodes_error("CVodeQuadReInitB", flag);

    } else {

      // Initialize the adjoint integration
      initAdj(m);
    }
    casadi_msg("CvodesInterface::resetB end");
  }

  void CvodesInterface::retreat(IntegratorMemory& mem, double t,
                                double* rx, double* rz, double* rq) const {
    CvodesMemory& m = dynamic_cast<CvodesMemory&>(mem);

    // Integrate, unless already at desired time
    if (t<m.t) {
      int flag = CVodeB(m.mem, t, CV_NORMAL);
      if (flag<CV_SUCCESS) cvodes_error("CVodeB", flag);

      // Get backward state
      flag = CVodeGetB(m.mem, m.whichB, &m.t, m.rxz);
      if (flag!=CV_SUCCESS) cvodes_error("CVodeGetB", flag);

      // Get backward qudratures
      if (nrq_>0) {
        flag = CVodeGetQuadB(m.mem, m.whichB, &m.t, m.rq);
        if (flag!=CV_SUCCESS) cvodes_error("CVodeGetQuadB", flag);
      }
    }

    // Save outputs
    casadi_copy(NV_DATA_S(m.rxz), nrx_, rx);
    casadi_copy(NV_DATA_S(m.rq), nrq_, rq);

    if (gather_stats_) {
      CVodeMem cv_mem = static_cast<CVodeMem>(m.mem);
      CVadjMem ca_mem = cv_mem->cv_adj_mem;
      CVodeBMem cvB_mem = ca_mem->cvB_mem;

      int flag = CVodeGetIntegratorStats(cvB_mem->cv_mem, &m.nstepsB,
                                         &m.nfevalsB, &m.nlinsetupsB, &m.netfailsB, &m.qlastB,
                                         &m.qcurB, &m.hinusedB, &m.hlastB, &m.hcurB, &m.tcurB);
      if (flag!=CV_SUCCESS) cvodes_error("CVodeGetIntegratorStats", flag);
    }
  }

  void CvodesInterface::printStats(IntegratorMemory& mem, std::ostream &stream) const {
    CvodesMemory& m = dynamic_cast<CvodesMemory&>(mem);

    long nsteps, nfevals, nlinsetups, netfails;
    int qlast, qcur;
    double hinused, hlast, hcur, tcur;
    int flag = CVodeGetIntegratorStats(m.mem, &nsteps, &nfevals, &nlinsetups, &netfails, &qlast,
                                       &qcur, &hinused, &hlast, &hcur, &tcur);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeGetIntegratorStats", flag);

    // Get the number of right hand side evaluations in the linear solver
    long nfevals_linsol=0;
    switch (linsol_f_) {
    case SD_DENSE:
    case SD_BANDED:
      flag = CVDlsGetNumRhsEvals(m.mem, &nfevals_linsol);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsGetNumRhsEvals", flag);
      break;
    case SD_ITERATIVE:
      flag = CVSpilsGetNumRhsEvals(m.mem, &nfevals_linsol);
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

    stream << "number of checkpoints stored: " << m.ncheck << endl;
    stream << std::endl;

    stream << "Time spent in the ODE residual: " << m.t_res << " s." << endl;
    stream << "Time spent in the forward sensitivity residual: " << m.t_fres << " s." << endl;
    stream << "Time spent in the jacobian function or jacobian times vector function: "
           << m.t_jac << " s." << endl;
    stream << "Time spent in the linear solver solve function: " << m.t_lsolve << " s." << endl;
    stream << "Time spent to generate the jacobian in the linear solver setup function: "
           << m.t_lsetup_jac << " s." << endl;
    stream << "Time spent to factorize the jacobian in the linear solver setup function: "
           << m.t_lsetup_fac << " s." << endl;
    stream << std::endl;
  }

  void CvodesInterface::cvodes_error(const string& module, int flag) {
    stringstream ss;
    ss << "Module \"" << module << "\" returned \"" << CVodeGetReturnFlagName(flag) << "\".";
    ss << " Consult Cvodes documentation.";
    casadi_error(ss.str());
  }

  void CvodesInterface::ehfun_wrapper(int error_code, const char *module, const char *function,
                                     char *msg, void *user_data) {
    try {
      casadi_assert(user_data);
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      m->self.ehfun(*m, error_code, module, function, msg);
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "ehfun failed: " << e.what() << endl;
    }
  }

  void CvodesInterface::ehfun(CvodesMemory& m, int error_code, const char *module,
                              const char *function, char *msg) const {
    if (!disable_internal_warnings_) {
      userOut<true, PL_WARN>() << msg << endl;
    }
  }

  void CvodesInterface::rhsS(CvodesMemory& m, int Ns, double t, N_Vector x, N_Vector xdot,
                             N_Vector *xF, N_Vector *xdotF, N_Vector tmp1, N_Vector tmp2) const {
    //    casadi_assert(Ns==nfdir_);

    // Record the current cpu time
    m.time1 = clock();

    // Commented out since a new implementation currently cannot be tested
    casadi_error("Commented out, #884, #794.");

    // Record timings
    m.time2 = clock();
    m.t_fres += static_cast<double>(m.time2-m.time1)/CLOCKS_PER_SEC;
  }

  int CvodesInterface::rhsS_wrapper(int Ns, double t, N_Vector x, N_Vector xdot, N_Vector *xF,
                                    N_Vector *xdotF, void *user_data,
                                    N_Vector tmp1, N_Vector tmp2) {
    try {
      casadi_assert(user_data);
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      m->self.rhsS(*m, Ns, t, x, xdot, xF, xdotF, tmp1, tmp2);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "fs failed: " << e.what() << endl;
      return 1;
    }
  }

  void CvodesInterface::rhsS1(CvodesMemory& m, int Ns, double t, N_Vector x, N_Vector xdot, int iS,
                              N_Vector xF, N_Vector xdotF, N_Vector tmp1, N_Vector tmp2) const {

    // Commented out since a new implementation currently cannot be tested
    casadi_error("Commented out, #884, #794.");
  }

  int CvodesInterface::rhsS1_wrapper(int Ns, double t, N_Vector x, N_Vector xdot, int iS,
                                    N_Vector xF, N_Vector xdotF, void *user_data,
                                    N_Vector tmp1, N_Vector tmp2) {
    try {
      casadi_assert(user_data);
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      m->self.rhsS1(*m, Ns, t, x, xdot, iS, xF, xdotF, tmp1, tmp2);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "fs failed: " << e.what() << endl;
      return 1;
    }
  }

  int CvodesInterface::rhsQ_wrapper(double t, N_Vector x, N_Vector qdot, void *user_data) {
    try {
      casadi_assert(user_data);
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      m->self.rhsQ(*m, t, x, qdot);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhsQ failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::rhsQ(CvodesMemory& m, double t, N_Vector x, N_Vector qdot) const {
    // Evaluate f_
    m.arg[DAE_T] = &t;
    m.arg[DAE_X] = NV_DATA_S(x);
    m.arg[DAE_Z] = 0;
    m.arg[DAE_P] = get_ptr(m.p);
    m.res[DAE_ODE] = 0;
    m.res[DAE_ALG] = 0;
    m.res[DAE_QUAD] = NV_DATA_S(qdot);
    f_(m.arg, m.res, m.iw, m.w, 0);
  }

  void CvodesInterface::rhsQS(CvodesMemory& m, int Ns, double t, N_Vector x, N_Vector *xF,
                              N_Vector qdot, N_Vector *qdotF, N_Vector tmp1, N_Vector tmp2) const {
    // Commented out since a new implementation currently cannot be tested
    casadi_error("Commented out, #884, #794.");
  }

  int CvodesInterface::rhsQS_wrapper(int Ns, double t, N_Vector x, N_Vector *xF, N_Vector qdot,
                                    N_Vector *qdotF, void *user_data,
                                    N_Vector tmp1, N_Vector tmp2) {
    try {
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      if (!m) {
        // SUNDIALS BUG!!!
        for (int i=0; i<Ns; ++i) N_VConst(0.0, qdotF[i]);
        return 0;
      }
      m->self.rhsQS(*m, Ns, t, x, xF, qdot, qdotF, tmp1, tmp2);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhsQS failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::
  rhsB(CvodesMemory& m, double t, N_Vector x, N_Vector rx, N_Vector rxdot) const {
    log("CvodesInterface::rhsB", "begin");

    // Debug output
    if (monitor_rhsB_) {
      printvar("t", t);
      printvar("x", x);
      printvar("rx", rx);
    }

    // Evaluate g_
    m.arg[RDAE_T] = &t;
    m.arg[RDAE_X] = NV_DATA_S(x);
    m.arg[RDAE_Z] = 0;
    m.arg[RDAE_P] = get_ptr(m.p);
    m.arg[RDAE_RX] = NV_DATA_S(rx);
    m.arg[RDAE_RZ] = 0;
    m.arg[RDAE_RP] = get_ptr(m.rp);
    m.res[RDAE_ODE] = NV_DATA_S(rxdot);
    m.res[RDAE_ALG] = 0;
    m.res[RDAE_QUAD] = 0;
    g_(m.arg, m.res, m.iw, m.w, 0);

    // Debug output
    if (monitor_rhsB_) {
      printvar("rxdot", rxdot);
    }

    // Negate (note definition of g)
    casadi_scal(nrx_, -1., NV_DATA_S(rxdot));

    log("CvodesInterface::rhsB", "end");
  }

  void CvodesInterface::rhsBS(CvodesMemory& m, double t, N_Vector x, N_Vector *xF, N_Vector rx,
                              N_Vector rxdot) const {

    // Commented out since a new implementation currently cannot be tested
    casadi_error("Commented out, #884, #794.");
  }

  int CvodesInterface::rhsB_wrapper(double t, N_Vector x, N_Vector rx, N_Vector rxdot,
                                   void *user_data) {
    try {
      casadi_assert(user_data);
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      m->self.rhsB(*m, t, x, rx, rxdot);
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
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      m->self.rhsBS(*m, t, x, xF, xB, xdotB);
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
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      m->self.rhsQB(*m, t, x, rx, rqdot);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhsQB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::rhsQB(CvodesMemory& m, double t, N_Vector x, N_Vector rx,
                              N_Vector rqdot) const {
    if (monitor_rhsQB_) {
      userOut() << "CvodesInterface::rhsQB: begin" << endl;
    }

    // Debug output
    if (monitor_rhsQB_) {
      printvar("t", t);
      printvar("x", x);
      printvar("rx", rx);
    }

    // Evaluate g_
    m.arg[RDAE_T] = &t;
    m.arg[RDAE_X] = NV_DATA_S(x);
    m.arg[RDAE_Z] = 0;
    m.arg[RDAE_P] = get_ptr(m.p);
    m.arg[RDAE_RX] = NV_DATA_S(rx);
    m.arg[RDAE_RZ] = 0;
    m.arg[RDAE_RP] = get_ptr(m.rp);
    m.res[RDAE_ODE] = 0;
    m.res[RDAE_ALG] = 0;
    m.res[RDAE_QUAD] = NV_DATA_S(rqdot);
    g_(m.arg, m.res, m.iw, m.w, 0);

    // Debug output
    if (monitor_rhsQB_) {
      printvar("rqdot", rqdot);
    }

    // Negate (note definition of g)
    casadi_scal(nrq_, -1., NV_DATA_S(rqdot));

    if (monitor_rhsQB_) {
      userOut() << "CvodesInterface::rhsQB: end" << endl;
    }
  }

  int CvodesInterface::jtimes_wrapper(N_Vector v, N_Vector Jv, double t, N_Vector x,
                                     N_Vector xdot, void *user_data, N_Vector tmp) {
    try {
      casadi_assert(user_data);
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      m->self.jtimes(*m, v, Jv, t, x, xdot, tmp);
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
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      m->self.jtimesB(*m, vB, JvB, t, x, xB, xdotB, tmpB);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "jtimes failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::jtimes(CvodesMemory& m, N_Vector v, N_Vector Jv, double t, N_Vector x,
                              N_Vector xdot, N_Vector tmp) const {
    log("CvodesInterface::jtimes", "begin");
    // Get time
    m.time1 = clock();

    // Evaluate f_fwd_
    m.arg[DAE_T] = &t;
    m.arg[DAE_X] = NV_DATA_S(x);
    m.arg[DAE_Z] = 0;
    m.arg[DAE_P] = get_ptr(m.p);
    m.arg[DAE_NUM_IN + DAE_T] = 0;
    m.arg[DAE_NUM_IN + DAE_X] = NV_DATA_S(v);
    m.arg[DAE_NUM_IN + DAE_Z] = 0;
    m.arg[DAE_NUM_IN + DAE_P] = 0;
    m.res[DAE_ODE] = 0;
    m.res[DAE_ALG] = 0;
    m.res[DAE_QUAD] = 0;
    m.res[DAE_NUM_OUT + DAE_ODE] = NV_DATA_S(Jv);
    m.res[DAE_NUM_OUT + DAE_ALG] = 0;
    m.res[DAE_NUM_OUT + DAE_QUAD] = 0;
    f_fwd_(m.arg, m.res, m.iw, m.w, 0);

    // Log time duration
    m.time2 = clock();
    m.t_jac += static_cast<double>(m.time2-m.time1)/CLOCKS_PER_SEC;

    log("CvodesInterface::jtimes", "end");
  }

  void CvodesInterface::jtimesB(CvodesMemory& m, N_Vector v, N_Vector Jv, double t, N_Vector x,
                                N_Vector rx, N_Vector rxdot, N_Vector tmpB) const {
    log("CvodesInterface::jtimesB", "begin");
    // Get time
    m.time1 = clock();

    // Hack:
    vector<const double*> arg1(g_fwd_.sz_arg());
    const double** arg1_ = get_ptr(arg1);
    vector<double*> res1(g_fwd_.sz_res());
    double** res1_ = get_ptr(res1);
    vector<int> iw(g_fwd_.sz_iw());
    int* iw_ = get_ptr(iw);
    vector<double> w(g_fwd_.sz_w());
    double* w_ = get_ptr(w);

    // Evaluate g_fwd_
    arg1_[RDAE_T] = &t;
    arg1_[RDAE_X] = NV_DATA_S(x);
    arg1_[RDAE_Z] = 0;
    arg1_[RDAE_P] = get_ptr(m.p);
    arg1_[RDAE_RX] = NV_DATA_S(rx);
    arg1_[RDAE_RZ] = 0;
    arg1_[RDAE_RP] = get_ptr(m.rp);;
    arg1_[RDAE_NUM_IN + RDAE_T] = 0;
    arg1_[RDAE_NUM_IN + RDAE_X] = 0;
    arg1_[RDAE_NUM_IN + RDAE_Z] = 0;
    arg1_[RDAE_NUM_IN + RDAE_P] = 0;
    arg1_[RDAE_NUM_IN + RDAE_RX] = NV_DATA_S(v);
    arg1_[RDAE_NUM_IN + RDAE_RZ] = 0;
    arg1_[RDAE_NUM_IN + RDAE_RP] = 0;
    res1_[RDAE_ODE] = 0;
    res1_[RDAE_ALG] = 0;
    res1_[RDAE_QUAD] = 0;
    res1_[RDAE_NUM_OUT + RDAE_ODE] = NV_DATA_S(Jv);
    res1_[RDAE_NUM_OUT + RDAE_ALG] = 0;
    res1_[RDAE_NUM_OUT + RDAE_QUAD] = 0;
    g_fwd_(arg1_, res1_, iw_, w_, 0);

    // Log time duration
    m.time2 = clock();
    m.t_jac += static_cast<double>(m.time2-m.time1)/CLOCKS_PER_SEC;
    log("CvodesInterface::jtimesB", "end");
  }

  int CvodesInterface::djac_wrapper(long N, double t, N_Vector x, N_Vector xdot, DlsMat Jac,
                                   void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    try {
      casadi_assert(user_data);
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      m->self.djac(*m, N, t, x, xdot, Jac, tmp1, tmp2, tmp3);
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
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      m->self.djacB(*m, NeqB, t, x, xB, xdotB, JacB, tmp1B, tmp2B, tmp3B);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "djacB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::djac(CvodesMemory& m, long N, double t, N_Vector x, N_Vector xdot,
                             DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) const {
    log("CvodesInterface::djac", "begin");

    // Get time
    m.time1 = clock();

    // Evaluate jac_
    m.arg[DAE_T] = &t;
    m.arg[DAE_X] = NV_DATA_S(x);
    m.arg[DAE_Z] = 0;
    m.arg[DAE_P] = get_ptr(m.p);
    double one=1, zero=0;
    m.arg[DAE_NUM_IN] = &one;
    m.arg[DAE_NUM_IN+1] = &zero;
    fill_n(m.res, jac_.n_out(), static_cast<double*>(0));
    m.res[0] = m.w + jac_.sz_w();
    jac_(m.arg, m.res, m.iw, m.w, 0);
    double *val = m.res[0];

    // Get sparsity and non-zero elements
    const int* colind = jac_.sparsity_out(0).colind();
    int ncol = jac_.size2_out(0);
    const int* row = jac_.sparsity_out(0).row();

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
    m.time2 = clock();
    m.t_jac += static_cast<double>(m.time2-m.time1)/CLOCKS_PER_SEC;

    log("CvodesInterface::djac", "end");
  }

  void CvodesInterface::djacB(CvodesMemory& m, long NeqB, double t, N_Vector x, N_Vector xB,
                              N_Vector xdotB, DlsMat JacB, N_Vector tmp1B, N_Vector tmp2B,
                              N_Vector tmp3B) const {
    log("CvodesInterface::djacB", "begin");
    // Get time
    m.time1 = clock();

    // Evaluate jacB_
    m.arg[RDAE_T] = &t;
    m.arg[RDAE_X] = NV_DATA_S(x);
    m.arg[RDAE_Z] = 0;
    m.arg[RDAE_P] = get_ptr(m.p);
    m.arg[RDAE_RX] = NV_DATA_S(xB);
    m.arg[RDAE_RZ] = 0;
    m.arg[RDAE_RP] = get_ptr(m.rp);
    double minus_one = -1;
    m.arg[RDAE_NUM_IN] = &minus_one;
    m.arg[RDAE_NUM_IN+1] = 0;
    fill_n(m.res, jacB_.n_out(), static_cast<double*>(0));
    m.res[0] = m.w + jacB_.sz_w();
    jacB_(m.arg, m.res, m.iw, m.w, 0);

    // Get sparsity and non-zero elements
    const int* colind = jacB_.sparsity_out(0).colind();
    int ncol = jacB_.size2_out(0);
    const int* row = jacB_.sparsity_out(0).row();
    double *val = m.res[0];

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
    m.time2 = clock();
    m.t_jac += static_cast<double>(m.time2-m.time1)/CLOCKS_PER_SEC;
    log("CvodesInterface::djacB", "end");
  }

  int CvodesInterface::bjac_wrapper(long N, long mupper, long mlower, double t, N_Vector x,
                                   N_Vector xdot, DlsMat Jac, void *user_data,
                                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    try {
      casadi_assert(user_data);
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      m->self.bjac(*m, N, mupper, mlower, t, x, xdot, Jac, tmp1, tmp2, tmp3);
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
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      m->self.bjacB(*m, NeqB, mupperB, mlowerB, t, x, xB, xdotB, JacB, tmp1B, tmp2B, tmp3B);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "bjacB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::bjac(CvodesMemory& m, long N, long mupper, long mlower, double t,
                             N_Vector x, N_Vector xdot,
                             DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) const {
    log("CvodesInterface::bjac", "begin");

    // Get time
    m.time1 = clock();

    // Evaluate jac_
    m.arg[DAE_T] = &t;
    m.arg[DAE_X] = NV_DATA_S(x);
    m.arg[DAE_Z] = 0;
    m.arg[DAE_P] = get_ptr(m.p);
    double one=1;
    m.arg[DAE_NUM_IN] = &one;
    m.arg[DAE_NUM_IN+1] = 0;
    fill_n(m.res, jac_.n_out(), static_cast<double*>(0));
    m.res[0] = m.w + jac_.sz_w();
    jac_(m.arg, m.res, m.iw, m.w, 0);
    double *val = m.res[0];

    // Get sparsity and non-zero elements
    const int* colind = jac_.sparsity_out(0).colind();
    int ncol = jac_.size2_out(0);
    const int* row = jac_.sparsity_out(0).row();

    // Loop over columns
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
    m.time2 = clock();
    m.t_jac += static_cast<double>(m.time2-m.time1)/CLOCKS_PER_SEC;

    log("CvodesInterface::bjac", "end");
  }

  void CvodesInterface::bjacB(CvodesMemory& m,
                              long NeqB, long mupperB, long mlowerB, double t, N_Vector x,
                              N_Vector xB, N_Vector xdotB, DlsMat JacB, N_Vector tmp1B,
                              N_Vector tmp2B, N_Vector tmp3B) const {
    log("CvodesInterface::bjacB", "begin");

    // Get time
    m.time1 = clock();

    // Evaluate jacB_
    m.arg[RDAE_T] = &t;
    m.arg[RDAE_X] = NV_DATA_S(x);
    m.arg[RDAE_Z] = 0;
    m.arg[RDAE_P] = get_ptr(m.p);
    m.arg[RDAE_RX] = NV_DATA_S(xB);
    m.arg[RDAE_RZ] = 0;
    m.arg[RDAE_RP] = get_ptr(m.rp);
    double minus_one = -1;
    m.arg[RDAE_NUM_IN] = &minus_one;
    m.arg[RDAE_NUM_IN+1] = 0;
    fill_n(m.res, jacB_.n_out(), static_cast<double*>(0));
    m.res[0] = m.w + jacB_.sz_w();
    jacB_(m.arg, m.res, m.iw, m.w, 0);

    // Get sparsity and non-zero elements
    const int* colind = jacB_.sparsity_out(0).colind();
    int ncol = jacB_.size2_out(0);
    const int* row = jacB_.sparsity_out(0).row();
    double *val = m.res[0];

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
    m.time2 = clock();
    m.t_jac += static_cast<double>(m.time2-m.time1)/CLOCKS_PER_SEC;

    log("CvodesInterface::bjacB", "end");
  }

  void CvodesInterface::setStopTime(IntegratorMemory& mem, double tf) const {
    // Set the stop time of the integration -- don't integrate past this point
    CvodesMemory& m = dynamic_cast<CvodesMemory&>(mem);
    int flag = CVodeSetStopTime(m.mem, tf);
    if (flag != CV_SUCCESS) cvodes_error("CVodeSetStopTime", flag);
  }

  int CvodesInterface::psolve_wrapper(double t, N_Vector x, N_Vector xdot, N_Vector r,
                                      N_Vector z, double gamma, double delta, int lr,
                                      void *user_data, N_Vector tmp) {
    try {
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      casadi_assert(m);
      m->self.psolve(*m, t, x, xdot, r, z, gamma, delta, lr, tmp);
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
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      casadi_assert(m);
      m->self.psolveB(*m, t, x, xB, xdotB, rvecB, zvecB, gammaB, deltaB, lr, tmpB);
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
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      casadi_assert(m);
      m->self.psetup(*m, t, x, xdot, jok, jcurPtr, gamma, tmp1, tmp2, tmp3);
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
      CvodesMemory* m = static_cast<CvodesMemory*>(user_data);
      casadi_assert(m);
      m->self.psetupB(*m, t, x, xB, xdotB, jokB, jcurPtrB, gammaB, tmp1B, tmp2B, tmp3B);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "psetupB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::psolve(CvodesMemory& m,
                               double t, N_Vector x, N_Vector xdot, N_Vector r, N_Vector z,
                               double gamma, double delta, int lr, N_Vector tmp) const {
    // Get time
    m.time1 = clock();

    // Copy input to output, if necessary
    if (r!=z) {
      N_VScale(1.0, r, z);
    }

    // Solve the (possibly factorized) system
    casadi_assert(linsol_.nnz_out(0) == NV_LENGTH_S(z));
    linsol_.linsol_solve(NV_DATA_S(z));

    // Log time duration
    m.time2 = clock();
    m.t_lsolve += static_cast<double>(m.time2-m.time1)/CLOCKS_PER_SEC;
  }

  void CvodesInterface::psolveB(CvodesMemory& m,
                                double t, N_Vector x, N_Vector xB, N_Vector xdotB, N_Vector rvecB,
                                N_Vector zvecB, double gammaB, double deltaB,
                                int lr, N_Vector tmpB) const {
    // Get time
    m.time1 = clock();

    // Copy input to output, if necessary
    if (rvecB!=zvecB) {
      N_VScale(1.0, rvecB, zvecB);
    }

    // Solve the (possibly factorized) system
    casadi_assert(linsolB_.nnz_out(0) == NV_LENGTH_S(zvecB));
    linsolB_.linsol_solve(NV_DATA_S(zvecB), 1);

    // Log time duration
    m.time2 = clock();
    m.t_lsolve += static_cast<double>(m.time2-m.time1)/CLOCKS_PER_SEC;
  }

  void CvodesInterface::psetup(CvodesMemory& m,
                               double t, N_Vector x, N_Vector xdot, booleantype jok,
                               booleantype *jcurPtr, double gamma, N_Vector tmp1,
                               N_Vector tmp2, N_Vector tmp3) const {
    log("CvodesInterface::psetup", "begin");
    // Get time
    m.time1 = clock();

    // Evaluate jac_
    m.arg[DAE_T] = &t;
    m.arg[DAE_X] = NV_DATA_S(x);
    m.arg[DAE_Z] = 0;
    m.arg[DAE_P] = get_ptr(m.p);
    double d1 = -gamma, d2 = 1.;
    m.arg[DAE_NUM_IN] = &d1;
    m.arg[DAE_NUM_IN+1] = &d2;
    fill_n(m.res, jac_.n_out(), static_cast<double*>(0));
    double *val = m.w;
    double *w2 = m.w + jac_.nnz_out(0);
    m.res[0] = val;
    jac_(m.arg, m.res, m.iw, w2, 0);

    // Log time duration
    m.time2 = clock();
    m.t_lsetup_jac += static_cast<double>(m.time2-m.time1)/CLOCKS_PER_SEC;

    // Prepare the solution of the linear system (e.g. factorize)
    linsol_.setup(m.arg+LINSOL_NUM_IN, m.res+LINSOL_NUM_OUT, m.iw, w2);
    linsol_.linsol_factorize(val);

    // Log time duration
    m.time1 = clock();
    m.t_lsetup_fac += static_cast<double>(m.time1-m.time2)/CLOCKS_PER_SEC;

    log("CvodesInterface::psetup", "end");
  }

  void CvodesInterface::psetupB(CvodesMemory& m,
                                double t, N_Vector x, N_Vector xB, N_Vector xdotB,
                                booleantype jokB, booleantype *jcurPtrB, double gammaB,
                                N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) const {
    log("CvodesInterface::psetupB", "begin");
    // Get time
    m.time1 = clock();

    // Evaluate jacB_
    m.arg[RDAE_T] = &t;
    m.arg[RDAE_X] = NV_DATA_S(x);
    m.arg[RDAE_Z] = 0;
    m.arg[RDAE_P] = get_ptr(m.p);
    m.arg[RDAE_RX] = NV_DATA_S(xB);
    m.arg[RDAE_RZ] = 0;
    m.arg[RDAE_RP] = get_ptr(m.rp);
    m.arg[RDAE_NUM_IN] = &gammaB;
    double one=1;
    m.arg[RDAE_NUM_IN+1] = &one;
    fill_n(m.res, jacB_.n_out(), static_cast<double*>(0));
    double *val = m.w;
    double *w2 = m.w + jacB_.nnz_out(0);
    m.res[0] = val;
    jacB_(m.arg, m.res, m.iw, w2, 0);

    // Log time duration
    m.time2 = clock();
    m.t_lsetup_jac += static_cast<double>(m.time2-m.time1)/CLOCKS_PER_SEC;

    // Prepare the solution of the linear system (e.g. factorize)
    linsolB_.setup(m.arg+LINSOL_NUM_IN, m.res+LINSOL_NUM_OUT, m.iw, w2);
    linsolB_.linsol_factorize(val);

    // Log time duration
    m.time1 = clock();
    m.t_lsetup_fac += static_cast<double>(m.time1-m.time2)/CLOCKS_PER_SEC;
    log("CvodesInterface::psetupB", "end");
  }

  void CvodesInterface::lsetup(CvodesMemory& m,
                               CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
                               booleantype *jcurPtr,
                               N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3) const {
    // Current time
    double t = cv_mem->cv_tn;

    // Scaling factor before J
    double gamma = cv_mem->cv_gamma;

    // Call the preconditioner setup function (which sets up the linear solver)
    psetup(m, t, x, xdot, FALSE, jcurPtr, gamma, vtemp1, vtemp2, vtemp3);
  }

  void CvodesInterface::lsetupB(CvodesMemory& m,
                                double t, double gamma, int convfail,
                                N_Vector x, N_Vector xB, N_Vector xdotB, booleantype *jcurPtr,
                                N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3) const {
    // Call the preconditioner setup function (which sets up the linear solver)
    psetupB(m, t, x, xB, xdotB, FALSE, jcurPtr, gamma, vtemp1, vtemp2, vtemp3);
  }

  int CvodesInterface::lsetup_wrapper(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
                                     booleantype *jcurPtr,
                                     N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3) {
    try {
      CvodesMemory* m = static_cast<CvodesMemory*>(cv_mem->cv_lmem);
      casadi_assert(m);
      m->self.lsetup(*m, cv_mem, convfail, x, xdot, jcurPtr, vtemp1, vtemp2, vtemp3);
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
      CvodesMemory* m = static_cast<CvodesMemory*>(cv_mem->cv_lmem);
      casadi_assert(m);
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

      m->self.lsetupB(*m, t, gamma, convfail, ca_mem->ca_ytmp, x, xdot,
                          jcurPtr, vtemp1, vtemp2, vtemp3);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "lsetupB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::lsolve(CvodesMemory& m,
                               CVodeMem cv_mem, N_Vector b, N_Vector weight,
                               N_Vector x, N_Vector xdot) const {
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
    psolve(m, t, x, xdot, b, b, gamma, delta, lr, 0);

    log("CvodesInterface::lsolve", "end");
  }

  void CvodesInterface::lsolveB(CvodesMemory& m,
                                double t, double gamma, N_Vector b, N_Vector weight,
                                N_Vector x, N_Vector xB, N_Vector xdotB) const {
    log("CvodesInterface::lsolveB", "begin");
    // Accuracy
    double delta = 0.0;

    // Left/right preconditioner
    int lr = 1;

    // Call the preconditioner solve function (which solves the linear system)
    psolveB(m, t, x, xB, xdotB, b, b, gamma, delta, lr, 0);

    log("CvodesInterface::lsolveB", "end");
  }

  int CvodesInterface::lsolve_wrapper(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                                     N_Vector x, N_Vector xdot) {
    try {
      CvodesMemory* m = static_cast<CvodesMemory*>(cv_mem->cv_lmem);
      casadi_assert(m);
      m->self.lsolve(*m, cv_mem, b, weight, x, xdot);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "lsolve failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::lsolveB_wrapper(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                                      N_Vector x, N_Vector xdot) {
    try {
      CvodesMemory* m = static_cast<CvodesMemory*>(cv_mem->cv_lmem);
      casadi_assert(m);
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

      m->self.lsolveB(*m, t, gamma, b, weight, ca_mem->ca_ytmp, x, xdot);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "lsolveB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::initDenseLinsol(CvodesMemory& m) const {
    int flag = CVDense(m.mem, nx_);
    if (flag!=CV_SUCCESS) cvodes_error("CVDense", flag);
    if (exact_jacobian_) {
      flag = CVDlsSetDenseJacFn(m.mem, djac_wrapper);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsSetDenseJacFn", flag);
    }
  }

  void CvodesInterface::initBandedLinsol(CvodesMemory& m) const {
    pair<int, int> bw = getBandwidth();
    int flag = CVBand(m.mem, nx_, bw.first, bw.second);
    if (flag!=CV_SUCCESS) cvodes_error("CVBand", flag);
    if (exact_jacobian_) {
      flag = CVDlsSetBandJacFn(m.mem, bjac_wrapper);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsSetBandJacFn", flag);
    }
  }

  void CvodesInterface::initIterativeLinsol(CvodesMemory& m) const {
    // Attach the sparse solver
    int flag;
    switch (itsol_f_) {
    case SD_GMRES:
      flag = CVSpgmr(m.mem, pretype_f_, max_krylov_);
      if (flag!=CV_SUCCESS) cvodes_error("CVBand", flag);
      break;
    case SD_BCGSTAB:
      flag = CVSpbcg(m.mem, pretype_f_, max_krylov_);
      if (flag!=CV_SUCCESS) cvodes_error("CVSpbcg", flag);
      break;
    case SD_TFQMR:
      flag = CVSptfqmr(m.mem, pretype_f_, max_krylov_);
      if (flag!=CV_SUCCESS) cvodes_error("CVSptfqmr", flag);
      break;
    }

    // Attach functions for jacobian information
    if (exact_jacobian_) {
      flag = CVSpilsSetJacTimesVecFn(m.mem, jtimes_wrapper);
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
      flag = CVSpilsSetPreconditioner(m.mem, psetup_wrapper, psolve_wrapper);
      if (flag != CV_SUCCESS) cvodes_error("CVSpilsSetPreconditioner", flag);
    }
  }

  void CvodesInterface::initUserDefinedLinsol(CvodesMemory& m) const {
    // Make sure that a Jacobian has been provided
    if (jac_.isNull())
        throw CasadiException("CvodesInterface::initUserDefinedLinsol(): "
                              "No Jacobian has been provided.");

    // Make sure that a linear solver has been provided
    if (linsol_.isNull())
        throw CasadiException("CvodesInterface::initUserDefinedLinsol(): "
                              "No user defined linear solver has been provided.");

    //  Set fields in the IDA memory
    CVodeMem cv_mem = static_cast<CVodeMem>(m.mem);
    cv_mem->cv_lmem   = &m;
    cv_mem->cv_lsetup = lsetup_wrapper;
    cv_mem->cv_lsolve = lsolve_wrapper;
    cv_mem->cv_setupNonNull = TRUE;
  }

  void CvodesInterface::initDenseLinsolB(CvodesMemory& m) const {
    int flag = CVDenseB(m.mem, m.whichB, nrx_);
    if (flag!=CV_SUCCESS) cvodes_error("CVDenseB", flag);
    if (exact_jacobianB_) {
      flag = CVDlsSetDenseJacFnB(m.mem, m.whichB, djacB_wrapper);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsSetDenseJacFnB", flag);
    }
  }

  void CvodesInterface::initBandedLinsolB(CvodesMemory& m) const {
    pair<int, int> bw = getBandwidthB();
    int flag = CVBandB(m.mem, m.whichB, nrx_, bw.first, bw.second);
    if (flag!=CV_SUCCESS) cvodes_error("CVBandB", flag);

    if (exact_jacobianB_) {
      flag = CVDlsSetBandJacFnB(m.mem, m.whichB, bjacB_wrapper);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsSetBandJacFnB", flag);
    }
  }

  void CvodesInterface::initIterativeLinsolB(CvodesMemory& m) const {
    int flag;
    switch (itsol_g_) {
    case SD_GMRES:
      flag = CVSpgmrB(m.mem, m.whichB, pretype_g_, max_krylovB_);
      if (flag!=CV_SUCCESS) cvodes_error("CVSpgmrB", flag);
      break;
    case SD_BCGSTAB:
      flag = CVSpbcgB(m.mem, m.whichB, pretype_g_, max_krylovB_);
      if (flag!=CV_SUCCESS) cvodes_error("CVSpbcgB", flag);
      break;
    case SD_TFQMR:
      flag = CVSptfqmrB(m.mem, m.whichB, pretype_g_, max_krylovB_);
      if (flag!=CV_SUCCESS) cvodes_error("CVSptfqmrB", flag);
      break;
    }

    // Attach functions for jacobian information
    if (exact_jacobianB_) {
      flag = CVSpilsSetJacTimesVecFnB(m.mem, m.whichB, jtimesB_wrapper);
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
      flag = CVSpilsSetPreconditionerB(m.mem, m.whichB, psetupB_wrapper, psolveB_wrapper);
      if (flag != CV_SUCCESS) cvodes_error("CVSpilsSetPreconditionerB", flag);
    }

  }

  void CvodesInterface::initUserDefinedLinsolB(CvodesMemory& m) const {
    // Make sure that a Jacobian has been provided
    if (jacB_.isNull())
        throw CasadiException("CvodesInterface::initUserDefinedLinsolB(): "
                              "No backwards Jacobian has been provided.");

    // Make sure that a linear solver has been provided
    if (linsolB_.isNull())
        throw CasadiException("CvodesInterface::initUserDefinedLinsolB(): "
                              "No user defined backward linear solver has been provided.");

    CVodeMem cv_mem = static_cast<CVodeMem>(m.mem);
    CVadjMem ca_mem = cv_mem->cv_adj_mem;
    CVodeBMem cvB_mem = ca_mem->cvB_mem;
    cvB_mem->cv_lmem   = &m;

    cvB_mem->cv_mem->cv_lmem = &m;
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

  CvodesMemory::CvodesMemory(const CvodesInterface& s) : self(s) {
    this->mem = 0;
    this->isInitAdj = false;

    // Reset checkpoints counter
    this->ncheck = 0;
  }

  CvodesMemory::~CvodesMemory() {
    if (this->mem) CVodeFree(&this->mem);
  }

  Dict CvodesMemory::get_stats() const {
    Dict stats = SundialsMemory::get_stats();
    return stats;
  }

} // namespace casadi
