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
    plugin->version = 30;
    return 0;
  }

  extern "C"
  void CASADI_INTEGRATOR_CVODES_EXPORT casadi_load_integrator_cvodes() {
    Integrator::registerPlugin(casadi_register_integrator_cvodes);
  }

  CvodesInterface::CvodesInterface(const std::string& name, const Function& dae)
    : SundialsInterface(name, dae) {
  }

  CvodesInterface::~CvodesInterface() {
    clear_memory();
  }

  Options CvodesInterface::options_
  = {{&SundialsInterface::options_},
     {{"linear_multistep_method",
       {OT_STRING,
        "Integrator scheme: BDF|adams"}},
      {"nonlinear_solver_iteration",
       {OT_STRING,
        "Nonlinear solver type: NEWTON|functional"}},
      {"fsens_all_at_once",
       {OT_BOOL,
        "Calculate all right hand sides of the sensitivity equations at once"}}
     }
  };

  void CvodesInterface::init(const Dict& opts) {
    log("CvodesInterface::init", "begin");

    // Initialize the base classes
    SundialsInterface::init(opts);

    // Default options
    string linear_multistep_method = "bdf";
    string nonlinear_solver_iteration = "newton";

    // Read options
    for (auto&& op : opts) {
      if (op.first=="linear_multistep_method") {
        linear_multistep_method = op.second.to_string();
      } else if (op.first=="nonlinear_solver_iteration") {
        nonlinear_solver_iteration = op.second.to_string();
      }
    }

    create_function("rhs", {"x", "p", "t"}, {"ode"});
    create_function("rhsQ", {"x", "p", "t"}, {"quad"});
    create_function("rhsB", {"rx", "rp", "x", "p", "t"},
                            {"rode"});
    create_function("rhsQB", {"rx", "rp", "x", "p", "t"},
                             {"rquad"});

    // Create a Jacobian if requested
    if (exact_jacobian_) {
      set_function(getJac());
      casadi_assert(get_function("jacF").sparsity_out(0) == sp_jac_dae_);
    }

    // Create a backwards Jacobian if requested
    if (exact_jacobianB_ && nrx_>0) {
      set_function(getJacB());
      casadi_assert(get_function("jacB").sparsity_out(0) == sp_jac_rdae_);
    };

    // Algebraic variables not supported
    casadi_assert_message(nz_==0 && nrz_==0,
                          "CVODES does not support algebraic variables");

    if (linear_multistep_method=="adams") {
      lmm_ = CV_ADAMS;
    } else if (linear_multistep_method=="bdf") {
      lmm_ = CV_BDF;
    } else {
      casadi_error("Unknown linear multistep method: " + linear_multistep_method);
    }

    if (nonlinear_solver_iteration=="newton") {
      iter_ = CV_NEWTON;
    } else if (nonlinear_solver_iteration=="functional") {
      iter_ = CV_FUNCTIONAL;
    } else {
      casadi_error("Unknown nonlinear solver iteration: " + nonlinear_solver_iteration);
    }

    // Attach functions for jacobian information
    if (exact_jacobian_) {
      switch (linsol_f_) {
      case SD_ITERATIVE:
        create_function("jtimes", {"t", "x", "p", "fwd:x"}, {"fwd:ode"});
        break;
      default: break;
      }
    }

    if (exact_jacobianB_) {
      switch (linsol_g_) {
      case SD_ITERATIVE:
        create_function("jtimesB",
          {"t", "x", "p", "rx", "rp", "fwd:rx"}, {"fwd:rode"});
        break;
      default: break;
      }
    }
  }

  void CvodesInterface::initAdj(CvodesMemory* m) const {

    // Create backward problem (use the same lmm and iter)
    int flag = CVodeCreateB(m->mem, lmm_, iter_, &m->whichB);
    if (flag != CV_SUCCESS) cvodes_error("CVodeCreateB", flag);

    // Initialize the backward problem
    double tB0 = grid_.back();
    flag = CVodeInitB(m->mem, m->whichB, rhsB_wrapper, tB0, m->rxz);
    if (flag != CV_SUCCESS) cvodes_error("CVodeInitB", flag);

    // Set tolerances
    flag = CVodeSStolerancesB(m->mem, m->whichB, reltolB_, abstolB_);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeSStolerancesB", flag);

    // User data
    flag = CVodeSetUserDataB(m->mem, m->whichB, m);
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
    flag = CVodeQuadInitB(m->mem, m->whichB, rhsQB_wrapper, m->rq);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeQuadInitB", flag);

    if (quad_err_con_) {
      flag = CVodeSetQuadErrConB(m->mem, m->whichB, true);
      if (flag != CV_SUCCESS) cvodes_error("CVodeSetQuadErrConB", flag);

      flag = CVodeQuadSStolerancesB(m->mem, m->whichB, reltolB_, abstolB_);
      if (flag != CV_SUCCESS) cvodes_error("CVodeQuadSStolerancesB", flag);
    }

    // Mark initialized
    m->isInitAdj = true;
  }

  void CvodesInterface::init_memory(void* mem) const {
    SundialsInterface::init_memory(mem);
    auto m = to_mem(mem);

    // Create CVodes memory block
    m->mem = CVodeCreate(lmm_, iter_);
    casadi_assert_message(m->mem!=0, "CVodeCreate: Creation failed");

    // Set error handler function
    int flag = CVodeSetErrHandlerFn(m->mem, ehfun_wrapper, &m);
    if (flag != CV_SUCCESS) cvodes_error("CVodeSetErrHandlerFn", flag);

    // Initialize CVodes
    double t0 = 0;
    flag = CVodeInit(m->mem, rhs_wrapper, t0, m->xz);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeInit", flag);

    // Set tolerances
    flag = CVodeSStolerances(m->mem, reltol_, abstol_);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeInit", flag);

    // Maximum number of steps
    CVodeSetMaxNumSteps(m->mem, max_num_steps_);
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
    flag = CVodeSetUserData(m->mem, m);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeSetUserData", flag);

    // Quadrature equations
    if (nq_>0) {
      // Initialize quadratures in CVodes
      flag = CVodeQuadInit(m->mem, rhsQ_wrapper, m->q);
      if (flag != CV_SUCCESS) cvodes_error("CVodeQuadInit", flag);

      // Should the quadrature errors be used for step size control?
      if (quad_err_con_) {
        flag = CVodeSetQuadErrCon(m->mem, true);
        if (flag != CV_SUCCESS) cvodes_error("CVodeSetQuadErrCon", flag);

        // Quadrature error tolerances
        // TODO(Joel): vector absolute tolerances
        flag = CVodeQuadSStolerances(m->mem, reltol_, abstol_);
        if (flag != CV_SUCCESS) cvodes_error("CVodeQuadSStolerances", flag);
      }
    }

    // Adjoint sensitivity problem
    if (nrx_>0) {
      // Get the interpolation type
      int interpType;
      if (interpolation_type_=="hermite") {
        interpType = CV_HERMITE;
      } else if (interpolation_type_=="polynomial") {
        interpType = CV_POLYNOMIAL;
      } else {
        casadi_error("\"interpolation_type\" must be \"hermite\" or \"polynomial\"");
      }

      // Initialize adjoint sensitivities
      flag = CVodeAdjInit(m->mem, steps_per_checkpoint_, interpType);
      if (flag != CV_SUCCESS) cvodes_error("CVodeAdjInit", flag);
      m->isInitAdj = false;
    }
  }

  int CvodesInterface::rhs_wrapper(double t, N_Vector x, N_Vector xdot, void *user_data) {
    try {
      casadi_assert(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      s.calc_function(m, "rhs", {NV_DATA_S(x), get_ptr(m->p), &t},
                              {NV_DATA_S(xdot)});
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhs failed: " << e.what() << endl;
      return 1;
    }
  }

  void CvodesInterface::reset(IntegratorMemory* mem, double t, const double* x,
                              const double* z, const double* _p) const {
    casadi_msg("CvodesInterface::reset begin");
    auto m = to_mem(mem);

    // Reset the base classes
    SundialsInterface::reset(mem, t, x, z, _p);

    // Re-initialize
    int flag = CVodeReInit(m->mem, t, m->xz);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeReInit", flag);

    // Re-initialize quadratures
    if (nq_>0) {
      N_VConst(0.0, m->q);
      flag = CVodeQuadReInit(m->mem, m->q);
      if (flag != CV_SUCCESS) cvodes_error("CVodeQuadReInit", flag);
    }

    // Turn off sensitivities
    flag = CVodeSensToggleOff(m->mem);
    if (flag != CV_SUCCESS) cvodes_error("CVodeSensToggleOff", flag);

    // Re-initialize backward integration
    if (nrx_>0) {
      flag = CVodeAdjReInit(m->mem);
      if (flag != CV_SUCCESS) cvodes_error("CVodeAdjReInit", flag);
    }

    // Set the stop time of the integration -- don't integrate past this point
    if (stop_at_end_) setStopTime(m, grid_.back());
    casadi_msg("CvodesInterface::reset end");
  }

  void CvodesInterface::advance(IntegratorMemory* mem, double t, double* x,
                                double* z, double* q) const {
    auto m = to_mem(mem);

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
    if (fabs(m->t-t)>=ttol) {
      // Integrate forward ...
      if (nrx_>0) {
        // ... with taping
        int flag = CVodeF(m->mem, t, m->xz, &m->t, CV_NORMAL, &m->ncheck);
        if (flag!=CV_SUCCESS && flag!=CV_TSTOP_RETURN) cvodes_error("CVodeF", flag);
      } else {
        // ... without taping
        int flag = CVode(m->mem, t, m->xz, &m->t, CV_NORMAL);
        if (flag!=CV_SUCCESS && flag!=CV_TSTOP_RETURN) cvodes_error("CVode", flag);
      }

      // Get quadratures
      if (nq_>0) {
        double tret;
        int flag = CVodeGetQuad(m->mem, &tret, m->q);
        if (flag!=CV_SUCCESS) cvodes_error("CVodeGetQuad", flag);
      }
    }

    // Set function outputs
    casadi_copy(NV_DATA_S(m->xz), nx_, x);
    casadi_copy(NV_DATA_S(m->q), nq_, q);

    // Print statistics
    if (print_stats_) printStats(m, userOut());

    int flag = CVodeGetIntegratorStats(m->mem, &m->nsteps, &m->nfevals, &m->nlinsetups,
                                       &m->netfails, &m->qlast, &m->qcur, &m->hinused,
                                       &m->hlast, &m->hcur, &m->tcur);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeGetIntegratorStats", flag);

    casadi_msg("CvodesInterface::integrate(" << t << ") end");
  }

  void CvodesInterface::resetB(IntegratorMemory* mem, double t, const double* rx,
                               const double* rz, const double* rp) const {
    auto m = to_mem(mem);

    // Reset the base classes
    SundialsInterface::resetB(mem, t, rx, rz, rp);

    int flag;
    if (m->isInitAdj) {
      flag = CVodeReInitB(m->mem, m->whichB, grid_.back(), m->rxz);
      if (flag != CV_SUCCESS) cvodes_error("CVodeReInitB", flag);

      flag = CVodeQuadReInitB(m->mem, m->whichB, m->rq);
      if (flag!=CV_SUCCESS) cvodes_error("CVodeQuadReInitB", flag);

    } else {

      // Initialize the adjoint integration
      initAdj(m);
    }
    casadi_msg("CvodesInterface::resetB end");
  }

  void CvodesInterface::retreat(IntegratorMemory* mem, double t,
                                double* rx, double* rz, double* rq) const {
    auto m = to_mem(mem);

    // Integrate, unless already at desired time
    if (t<m->t) {
      int flag = CVodeB(m->mem, t, CV_NORMAL);
      if (flag<CV_SUCCESS) cvodes_error("CVodeB", flag);

      // Get backward state
      flag = CVodeGetB(m->mem, m->whichB, &m->t, m->rxz);
      if (flag!=CV_SUCCESS) cvodes_error("CVodeGetB", flag);

      // Get backward qudratures
      if (nrq_>0) {
        flag = CVodeGetQuadB(m->mem, m->whichB, &m->t, m->rq);
        if (flag!=CV_SUCCESS) cvodes_error("CVodeGetQuadB", flag);
      }
    }

    // Save outputs
    casadi_copy(NV_DATA_S(m->rxz), nrx_, rx);
    casadi_copy(NV_DATA_S(m->rq), nrq_, rq);

    CVodeMem cv_mem = static_cast<CVodeMem>(m->mem);
    CVadjMem ca_mem = cv_mem->cv_adj_mem;
    CVodeBMem cvB_mem = ca_mem->cvB_mem;

    int flag = CVodeGetIntegratorStats(cvB_mem->cv_mem, &m->nstepsB,
                                       &m->nfevalsB, &m->nlinsetupsB, &m->netfailsB, &m->qlastB,
                                       &m->qcurB, &m->hinusedB, &m->hlastB, &m->hcurB, &m->tcurB);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeGetIntegratorStats", flag);
  }

  void CvodesInterface::printStats(IntegratorMemory* mem, std::ostream &stream) const {
    auto m = to_mem(mem);

    long nsteps, nfevals, nlinsetups, netfails;
    int qlast, qcur;
    double hinused, hlast, hcur, tcur;
    int flag = CVodeGetIntegratorStats(m->mem, &nsteps, &nfevals, &nlinsetups, &netfails, &qlast,
                                       &qcur, &hinused, &hlast, &hcur, &tcur);
    if (flag!=CV_SUCCESS) cvodes_error("CVodeGetIntegratorStats", flag);

    // Get the number of right hand side evaluations in the linear solver
    long nfevals_linsol=0;
    switch (linsol_f_) {
    case SD_DENSE:
    case SD_BANDED:
      flag = CVDlsGetNumRhsEvals(m->mem, &nfevals_linsol);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsGetNumRhsEvals", flag);
      break;
    case SD_ITERATIVE:
      flag = CVSpilsGetNumRhsEvals(m->mem, &nfevals_linsol);
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

    stream << "number of checkpoints stored: " << m->ncheck << endl;
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
      auto m = to_mem(user_data);
      auto& s = m->self;
      if (!s.disable_internal_warnings_) {
        userOut<true, PL_WARN>() << msg << endl;
      }
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "ehfun failed: " << e.what() << endl;
    }
  }

  int CvodesInterface::rhsS_wrapper(int Ns, double t, N_Vector x, N_Vector xdot, N_Vector *xF,
                                    N_Vector *xdotF, void *user_data,
                                    N_Vector tmp1, N_Vector tmp2) {
    try {
      casadi_assert(user_data);
      auto m = to_mem(user_data);
      casadi_error("Commented out, #884, #794.");
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "fs failed: " << e.what() << endl;
      return 1;
    }
  }

  int CvodesInterface::rhsS1_wrapper(int Ns, double t, N_Vector x, N_Vector xdot, int iS,
                                    N_Vector xF, N_Vector xdotF, void *user_data,
                                    N_Vector tmp1, N_Vector tmp2) {
    try {
      casadi_assert(user_data);
      auto m = to_mem(user_data);
      casadi_error("Commented out, #884, #794.");
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "fs failed: " << e.what() << endl;
      return 1;
    }
  }

  int CvodesInterface::rhsQ_wrapper(double t, N_Vector x, N_Vector qdot, void *user_data) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;
      s.calc_function(m, "rhsQ", {NV_DATA_S(x), get_ptr(m->p), &t},
                               {NV_DATA_S(qdot)});
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhsQ failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::rhsQS_wrapper(int Ns, double t, N_Vector x, N_Vector *xF, N_Vector qdot,
                                    N_Vector *qdotF, void *user_data,
                                    N_Vector tmp1, N_Vector tmp2) {
    try {
      if (!user_data) {
        // SUNDIALS BUG!!!
        for (int i=0; i<Ns; ++i) N_VConst(0.0, qdotF[i]);
        return 0;
      }
      auto m = to_mem(user_data);
      casadi_error("Commented out, #884, #794.");
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhsQS failed: " << e.what() << endl;;
      return 1;
    }
  }


  int CvodesInterface::rhsB_wrapper(double t, N_Vector x, N_Vector rx, N_Vector rxdot,
                                   void *user_data) {
    try {
      casadi_assert(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      s.calc_function(m, "rhsB", {NV_DATA_S(rx), get_ptr(m->rp),
                                  NV_DATA_S(x), get_ptr(m->p), &t},
                                 {NV_DATA_S(rxdot)});

      // Negate (note definition of g)
      casadi_scal(s.nrx_, -1., NV_DATA_S(rxdot));

      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhsB failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::rhsBS_wrapper(double t, N_Vector x, N_Vector *xF, N_Vector xB,
                                    N_Vector xdotB, void *user_data) {
    try {
      auto m = to_mem(user_data);
      casadi_error("Commented out, #884, #794.");
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
      auto m = to_mem(user_data);
      auto& s = m->self;
      s.calc_function(m, "rhsQB", {NV_DATA_S(rx), get_ptr(m->rp),
                                 NV_DATA_S(x), get_ptr(m->p), &t},
                                {NV_DATA_S(rqdot)});

      // Negate (note definition of g)
      casadi_scal(s.nrq_, -1., NV_DATA_S(rqdot));

      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhsQB failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::jtimes_wrapper(N_Vector v, N_Vector Jv, double t, N_Vector x,
                                     N_Vector xdot, void *user_data, N_Vector tmp) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;
      s.calc_function(m, "jtimes", {&t, NV_DATA_S(x), get_ptr(m->p), NV_DATA_S(v)},
        {NV_DATA_S(Jv)});
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "jtimes failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::jtimesB_wrapper(N_Vector v, N_Vector Jv, double t, N_Vector x,
                                      N_Vector rx, N_Vector rxdot, void *user_data ,
                                      N_Vector tmpB) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;
      s.calc_function(m, "jtimesB",
        {&t, NV_DATA_S(x), get_ptr(m->p), NV_DATA_S(rx), get_ptr(m->rp), NV_DATA_S(v)},
        {NV_DATA_S(Jv)});
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "jtimes failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::djac_wrapper(long N, double t, N_Vector x, N_Vector xdot, DlsMat Jac,
                                   void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    try {
      casadi_assert(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      double one=1;
      s.calc_function(m, "jacF",
        {&t, NV_DATA_S(x), get_ptr(m->p), &one, 0},
        {m->jac});

      // Save to Jac
      const int* colind = s.sp_jac_dae_.colind();
      int ncol = s.sp_jac_dae_.size2();
      const int* row = s.sp_jac_dae_.row();
      for (int cc=0; cc<ncol; ++cc) {
        for (int el=colind[cc]; el<colind[cc+1]; ++el) {
          DENSE_ELEM(Jac, row[el], cc) = m->jac[el];
        }
      }
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "djac failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::djacB_wrapper(long NeqB, double t, N_Vector x, N_Vector rx, N_Vector rxdot,
                                    DlsMat JacB, void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                                    N_Vector tmp3B) {
    try {
      casadi_assert(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      double minus_one = -1;
      s.calc_function(m, "jacB", {&t, NV_DATA_S(rx), get_ptr(m->rp),
                                 NV_DATA_S(x), get_ptr(m->p), &minus_one, 0},
                                {m->jacB});

      // Save to JacB
      const int* colind = s.sp_jac_rdae_.colind();
      int ncol = s.sp_jac_rdae_.size2();
      const int* row = s.sp_jac_rdae_.row();
      for (int cc=0; cc<ncol; ++cc) {
        for (int el=colind[cc]; el<colind[cc+1]; ++el) {
          DENSE_ELEM(JacB, row[el], cc) = m->jacB[el];
        }
      }
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "djacB failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::bjac_wrapper(long N, long mupper, long mlower, double t, N_Vector x,
                                   N_Vector xdot, DlsMat Jac, void *user_data,
                                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    try {
      casadi_assert(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      double one=1;
      s.calc_function(m, "jacF", {&t, NV_DATA_S(x), get_ptr(m->p), &one, 0},
                              {m->jac});

      // Save to Jac
      const int* colind = s.sp_jac_dae_.colind();
      int ncol = s.sp_jac_dae_.size2();
      const int* row = s.sp_jac_dae_.row();
      for (int cc=0; cc<ncol; ++cc) {
        for (int el=colind[cc]; el<colind[cc+1]; ++el) {
          int rr = row[el];
          if (cc-rr<=mupper && rr-cc<=mlower)
            BAND_ELEM(Jac, rr, cc) = m->jac[el];
        }
      }
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "bjac failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::bjacB_wrapper(long NeqB, long mupperB, long mlowerB, double t, N_Vector x,
                                    N_Vector rx, N_Vector rxdot, DlsMat JacB, void *user_data,
                                    N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
    try {
      casadi_assert(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      double minus_one = -1;
      s.calc_function(m, "jacB", {&t, NV_DATA_S(rx), get_ptr(m->rp),
                                 NV_DATA_S(x), get_ptr(m->p), &minus_one, 0},
                                {m->jacB});

      // Save to JacB
      const int* colind = s.sp_jac_rdae_.colind();
      int ncol = s.sp_jac_rdae_.size2();
      const int* row = s.sp_jac_rdae_.row();
      for (int cc=0; cc<ncol; ++cc) {
        for (int el=colind[cc]; el<colind[cc+1]; ++el) {
          int rr = row[el];
          if (cc-rr<=mupperB && rr-cc<=mlowerB)
            BAND_ELEM(JacB, rr, cc) = m->jacB[el];
        }
      }
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "bjacB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::setStopTime(IntegratorMemory* mem, double tf) const {
    // Set the stop time of the integration -- don't integrate past this point
    auto m = to_mem(mem);
    int flag = CVodeSetStopTime(m->mem, tf);
    if (flag != CV_SUCCESS) cvodes_error("CVodeSetStopTime", flag);
  }

  int CvodesInterface::psolve_wrapper(double t, N_Vector x, N_Vector xdot, N_Vector r,
                                      N_Vector z, double gamma, double delta, int lr,
                                      void *user_data, N_Vector tmp) {
    try {
      auto m = to_mem(user_data);
      m->self.psolve(m, t, x, xdot, r, z, gamma, delta, lr, tmp);
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
      auto m = to_mem(user_data);
      m->self.psolveB(m, t, x, xB, xdotB, rvecB, zvecB, gammaB, deltaB, lr, tmpB);
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
      auto m = to_mem(user_data);
      m->self.psetup(m, t, x, xdot, jok, jcurPtr, gamma, tmp1, tmp2, tmp3);
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
      auto m = to_mem(user_data);
      m->self.psetupB(m, t, x, xB, xdotB, jokB, jcurPtrB, gammaB, tmp1B, tmp2B, tmp3B);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "psetupB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::psolve(CvodesMemory* m,
                               double t, N_Vector x, N_Vector xdot, N_Vector r, N_Vector z,
                               double gamma, double delta, int lr, N_Vector tmp) const {
    // Copy input to output, if necessary
    if (r!=z) {
      N_VScale(1.0, r, z);
    }

    // Solve the (possibly factorized) system
    casadi_assert(linsol_.nnz_out(0) == NV_LENGTH_S(z));
    linsol_.linsol_solve(NV_DATA_S(z));
  }

  void CvodesInterface::psolveB(CvodesMemory* m,
                                double t, N_Vector x, N_Vector xB, N_Vector xdotB, N_Vector rvecB,
                                N_Vector zvecB, double gammaB, double deltaB,
                                int lr, N_Vector tmpB) const {
    // Copy input to output, if necessary
    if (rvecB!=zvecB) {
      N_VScale(1.0, rvecB, zvecB);
    }

    // Solve the (possibly factorized) system
    casadi_assert(linsolB_.nnz_out(0) == NV_LENGTH_S(zvecB));
    linsolB_.linsol_solve(NV_DATA_S(zvecB), 1);
  }

  void CvodesInterface::psetup(CvodesMemory* m,
                               double t, N_Vector x, N_Vector xdot, booleantype jok,
                               booleantype *jcurPtr, double gamma, N_Vector tmp1,
                               N_Vector tmp2, N_Vector tmp3) const {
    log("CvodesInterface::psetup", "begin");
    // Calculate Jacobian
    double d1 = -gamma, d2 = 1.;
    calc_function(m, "jacF", {&t, NV_DATA_S(x), get_ptr(m->p), &d1, &d2},
                            {m->jac});

    // Prepare the solution of the linear system (e.g. factorize)
    linsol_.setup(m->arg+LINSOL_NUM_IN, m->res+LINSOL_NUM_OUT, m->iw, m->w);
    linsol_.linsol_factorize(m->jac);

    log("CvodesInterface::psetup", "end");
  }

  void CvodesInterface::psetupB(CvodesMemory* m,
                                double t, N_Vector x, N_Vector rx, N_Vector rxdot,
                                booleantype jokB, booleantype *jcurPtrB, double gammaB,
                                N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) const {
    log("CvodesInterface::psetupB", "begin");

    // Calculate Jacobian
    double one=1;
    calc_function(m, "jacB", {&t, NV_DATA_S(rx), get_ptr(m->rp),
                               NV_DATA_S(x), get_ptr(m->p), &gammaB, &one},
                              {m->jacB});

    // Prepare the solution of the linear system (e.g. factorize)
    linsolB_.setup(m->arg+LINSOL_NUM_IN, m->res+LINSOL_NUM_OUT, m->iw, m->w);
    linsolB_.linsol_factorize(m->jacB);

    log("CvodesInterface::psetupB", "end");
  }

  int CvodesInterface::lsetup_wrapper(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
                                     booleantype *jcurPtr,
                                     N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3) {
    try {
      auto m = to_mem(cv_mem->cv_lmem);
      auto& s = m->self;

      // Current time
      double t = cv_mem->cv_tn;

      // Scaling factor before J
      double gamma = cv_mem->cv_gamma;

      // Call the preconditioner setup function (which sets up the linear solver)
      s.psetup(m, t, x, xdot, FALSE, jcurPtr, gamma, vtemp1, vtemp2, vtemp3);

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
      auto m = to_mem(cv_mem->cv_lmem);
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

      // Call the preconditioner setup function (which sets up the linear solver)
      if (psetupB_wrapper(t, ca_mem->ca_ytmp, x, xdot, FALSE, jcurPtr,
        gamma, static_cast<void*>(m), vtemp1, vtemp2, vtemp3)) return 1;

      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "lsetupB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::lsolve(CvodesMemory* m,
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

  void CvodesInterface::lsolveB(CvodesMemory* m,
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
      auto m = to_mem(cv_mem->cv_lmem);
      m->self.lsolve(m, cv_mem, b, weight, x, xdot);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "lsolve failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::lsolveB_wrapper(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                                      N_Vector x, N_Vector xdot) {
    try {
      auto m = to_mem(cv_mem->cv_lmem);
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

      m->self.lsolveB(m, t, gamma, b, weight, ca_mem->ca_ytmp, x, xdot);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "lsolveB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::initDenseLinsol(CvodesMemory* m) const {
    int flag = CVDense(m->mem, nx_);
    if (flag!=CV_SUCCESS) cvodes_error("CVDense", flag);
    if (exact_jacobian_) {
      flag = CVDlsSetDenseJacFn(m->mem, djac_wrapper);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsSetDenseJacFn", flag);
    }
  }

  void CvodesInterface::initBandedLinsol(CvodesMemory* m) const {
    pair<int, int> bw = getBandwidth();
    int flag = CVBand(m->mem, nx_, bw.first, bw.second);
    if (flag!=CV_SUCCESS) cvodes_error("CVBand", flag);
    if (exact_jacobian_) {
      flag = CVDlsSetBandJacFn(m->mem, bjac_wrapper);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsSetBandJacFn", flag);
    }
  }

  void CvodesInterface::initIterativeLinsol(CvodesMemory* m) const {
    // Attach the sparse solver
    int flag;
    switch (itsol_f_) {
    case SD_GMRES:
      flag = CVSpgmr(m->mem, pretype_f_, max_krylov_);
      if (flag!=CV_SUCCESS) cvodes_error("CVBand", flag);
      break;
    case SD_BCGSTAB:
      flag = CVSpbcg(m->mem, pretype_f_, max_krylov_);
      if (flag!=CV_SUCCESS) cvodes_error("CVSpbcg", flag);
      break;
    case SD_TFQMR:
      flag = CVSptfqmr(m->mem, pretype_f_, max_krylov_);
      if (flag!=CV_SUCCESS) cvodes_error("CVSptfqmr", flag);
      break;
    }

    // Attach functions for jacobian information
    if (exact_jacobian_) {
      flag = CVSpilsSetJacTimesVecFn(m->mem, jtimes_wrapper);
      if (flag!=CV_SUCCESS) cvodes_error("CVSpilsSetJacTimesVecFn", flag);
    }

    // Add a preconditioner
    if (use_preconditioner_) {
      casadi_assert_message(has_function("jacF"), "No Jacobian function");

      // Make sure that a linear solver has been provided
      if (linsol_.is_null())
          throw CasadiException("CvodesInterface::init(): "
                                "No user defined linear solver has been provided.");

      // Pass to IDA
      flag = CVSpilsSetPreconditioner(m->mem, psetup_wrapper, psolve_wrapper);
      if (flag != CV_SUCCESS) cvodes_error("CVSpilsSetPreconditioner", flag);
    }
  }

  void CvodesInterface::initUserDefinedLinsol(CvodesMemory* m) const {
    casadi_assert_message(has_function("jacF"), "No Jacobian function");

    // Make sure that a linear solver has been provided
    if (linsol_.is_null())
        throw CasadiException("CvodesInterface::initUserDefinedLinsol(): "
                              "No user defined linear solver has been provided.");

    //  Set fields in the IDA memory
    CVodeMem cv_mem = static_cast<CVodeMem>(m->mem);
    cv_mem->cv_lmem   = m;
    cv_mem->cv_lsetup = lsetup_wrapper;
    cv_mem->cv_lsolve = lsolve_wrapper;
    cv_mem->cv_setupNonNull = TRUE;
  }

  void CvodesInterface::initDenseLinsolB(CvodesMemory* m) const {
    int flag = CVDenseB(m->mem, m->whichB, nrx_);
    if (flag!=CV_SUCCESS) cvodes_error("CVDenseB", flag);
    if (exact_jacobianB_) {
      flag = CVDlsSetDenseJacFnB(m->mem, m->whichB, djacB_wrapper);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsSetDenseJacFnB", flag);
    }
  }

  void CvodesInterface::initBandedLinsolB(CvodesMemory* m) const {
    pair<int, int> bw = getBandwidthB();
    int flag = CVBandB(m->mem, m->whichB, nrx_, bw.first, bw.second);
    if (flag!=CV_SUCCESS) cvodes_error("CVBandB", flag);

    if (exact_jacobianB_) {
      flag = CVDlsSetBandJacFnB(m->mem, m->whichB, bjacB_wrapper);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsSetBandJacFnB", flag);
    }
  }

  void CvodesInterface::initIterativeLinsolB(CvodesMemory* m) const {
    int flag;
    switch (itsol_g_) {
    case SD_GMRES:
      flag = CVSpgmrB(m->mem, m->whichB, pretype_g_, max_krylovB_);
      if (flag!=CV_SUCCESS) cvodes_error("CVSpgmrB", flag);
      break;
    case SD_BCGSTAB:
      flag = CVSpbcgB(m->mem, m->whichB, pretype_g_, max_krylovB_);
      if (flag!=CV_SUCCESS) cvodes_error("CVSpbcgB", flag);
      break;
    case SD_TFQMR:
      flag = CVSptfqmrB(m->mem, m->whichB, pretype_g_, max_krylovB_);
      if (flag!=CV_SUCCESS) cvodes_error("CVSptfqmrB", flag);
      break;
    }

    // Attach functions for jacobian information
    if (exact_jacobianB_) {
      flag = CVSpilsSetJacTimesVecFnB(m->mem, m->whichB, jtimesB_wrapper);
      if (flag!=CV_SUCCESS) cvodes_error("CVSpilsSetJacTimesVecFnB", flag);
    }

    // Add a preconditioner
    if (use_preconditionerB_) {
      casadi_assert_message(has_function("jacB"), "No Jacobian function");

      // Make sure that a linear solver has been provided
      if (linsolB_.is_null())
        casadi_error("CvodesInterface::init(): "
                     "No user defined backwards  linear solver has been provided.");

      // Pass to IDA
      flag = CVSpilsSetPreconditionerB(m->mem, m->whichB, psetupB_wrapper, psolveB_wrapper);
      if (flag != CV_SUCCESS) cvodes_error("CVSpilsSetPreconditionerB", flag);
    }

  }

  void CvodesInterface::initUserDefinedLinsolB(CvodesMemory* m) const {
    casadi_assert_message(has_function("jacB"), "No Jacobian function");

    // Make sure that a linear solver has been provided
    if (linsolB_.is_null())
        throw CasadiException("CvodesInterface::initUserDefinedLinsolB(): "
                              "No user defined backward linear solver has been provided.");

    CVodeMem cv_mem = static_cast<CVodeMem>(m->mem);
    CVadjMem ca_mem = cv_mem->cv_adj_mem;
    CVodeBMem cvB_mem = ca_mem->cvB_mem;
    cvB_mem->cv_lmem   = m;

    cvB_mem->cv_mem->cv_lmem = m;
    cvB_mem->cv_mem->cv_lsetup = lsetupB_wrapper;
    cvB_mem->cv_mem->cv_lsolve = lsolveB_wrapper;
    cvB_mem->cv_mem->cv_setupNonNull = TRUE;
  }

  template<typename MatType>
  Function CvodesInterface::getJacGen() {
    vector<MatType> a = MatType::get_input(oracle_);
    vector<MatType> r = oracle_(a);

    // Get the Jacobian in the Newton iteration
    MatType c_x = MatType::sym("c_x");
    MatType c_xdot = MatType::sym("c_xdot");
    MatType jac = c_x*MatType::jacobian(r[DE_ODE], a[DE_X]) + c_xdot*MatType::eye(nx_);

    // Return generated function
    return Function("jacF", {a[DE_T], a[DE_X], a[DE_P], c_x, c_xdot}, {jac});
  }

  template<typename MatType>
  Function CvodesInterface::getJacGenB() {
    vector<MatType> a = MatType::get_input(oracle_);
    vector<MatType> r = oracle_(a);

    // Get the Jacobian in the Newton iteration
    MatType c_x = MatType::sym("c_x");
    MatType c_xdot = MatType::sym("c_xdot");
    MatType jac = c_x*MatType::jacobian(r[DE_RODE], a[DE_RX]) + c_xdot*MatType::eye(nrx_);

    // return generated function
    return Function("jacB",
      {a[DE_T], a[DE_RX], a[DE_RP], a[DE_X], a[DE_P], c_x, c_xdot}, {jac});
  }

  Function CvodesInterface::getJacB() {
    if (oracle_.is_a("sxfunction")) {
      return getJacGenB<SX>();
    } else {
      return getJacGenB<MX>();
    }
  }


  Function CvodesInterface::getJac() {
    if (oracle_.is_a("sxfunction")) {
      return getJacGen<SX>();
    } else {
      return getJacGen<MX>();
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

  Dict CvodesInterface::get_stats(void* mem) const {
    Dict stats = SundialsInterface::get_stats(mem);
    return stats;
  }

} // namespace casadi
