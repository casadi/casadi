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

    if (ns_==0) {
      // Create functions
      create_function("odeF", {"t", "x", "p"}, {"ode"});
      create_function("quadF", {"t", "x", "p"}, {"quad"});
      create_function("odeB", {"t", "rx", "rp", "x", "p"}, {"rode"});
      create_function("quadB", {"t", "rx", "rp", "x", "p"}, {"rquad"});
    } else {
      // Fetch nondifferentiated functions
      Function f = derivative_of_.get_function("odeF");
      set_function(f, "odeF");
      set_function(f.forward(ns_), "odeFS");
      if (nq1_>0) {
        f = derivative_of_.get_function("quadF");
        set_function(f, "quadF");
        set_function(f.forward(ns_), "quadFS");
      }
      if (nrx1_>0) {
        f = derivative_of_.get_function("odeB");
        set_function(f, "odeB");
        set_function(f.forward(ns_), "odeBS");
      }
      if (nrq1_>0) {
        f = derivative_of_.get_function("quadB");
        set_function(f, "quadB");
        set_function(f.forward(ns_), "quadBS");
      }
    }

    // Create a Jacobian if requested
    if (exact_jacobian_) {
      if (ns_==0) {
        set_function(oracle_.is_a("sxfunction") ? getJacF<SX>() : getJacF<MX>(), "jacF");
      } else {
        set_function(derivative_of_.get_function("jacF"), "jacF");
      }
      init_linsol();
    }

    // Create a backwards Jacobian if requested
    if (exact_jacobianB_ && nrx1_>0) {
      if (ns_==0) {
        set_function(oracle_.is_a("sxfunction") ? getJacB<SX>() : getJacB<MX>(), "jacB");
      } else {
        set_function(derivative_of_.get_function("jacB"), "jacB");
      }
      init_linsolB();
    }

    // Algebraic variables not supported
    casadi_assert_message(nz1_==0 && nrz_==0,
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
        create_function("jtimesF", {"t", "x", "p", "fwd:x"}, {"fwd:ode"});
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
    flag = CVodeInitB(m->mem, m->whichB, rhsB, tB0, m->rxz);
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
    flag = CVodeQuadInitB(m->mem, m->whichB, rhsQB, m->rq);
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
    int flag = CVodeSetErrHandlerFn(m->mem, ehfun, &m);
    if (flag != CV_SUCCESS) cvodes_error("CVodeSetErrHandlerFn", flag);

    // Initialize CVodes
    double t0 = 0;
    flag = CVodeInit(m->mem, rhs, t0, m->xz);
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
    if (nq1_>0) {
      // Initialize quadratures in CVodes
      flag = CVodeQuadInit(m->mem, rhsQ, m->q);
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

    // Forward sensitivity analysis
    if (ns_>0) {
      // Initialize forward sensitivity analysis
      flag = CVodeSensInit(m->mem, ns_, CV_SIMULTANEOUS, rhsS, m->xzS);
      if (flag!=CV_SUCCESS) cvodes_error("CVodeSensInit", flag);

      // Set sensitivity tolerances
      vector<double> abstolv(ns_, abstol_);
      flag = CVodeSensSStolerances(m->mem, reltol_, get_ptr(abstolv));
      if (flag!=CV_SUCCESS) cvodes_error("CVodeSensSStolerances", flag);

      // Quadrature equations
      if (nq1_>0) {
        // Initialize quadratures in CVodes
        flag = CVodeQuadSensInit(m->mem, rhsQS, m->qS);
        if (flag != CV_SUCCESS) cvodes_error("CVodeQuadSensInit", flag);

        // Should the quadrature errors be used for step size control?
        if (quad_err_con_) {
          flag = CVodeSetQuadSensErrCon(m->mem, true);
          if (flag != CV_SUCCESS) cvodes_error("CVodeSetQuadErrCon", flag);

          // Quadrature error tolerances
          flag = CVodeQuadSensSStolerances(m->mem, reltol_, get_ptr(abstolv));
          if (flag != CV_SUCCESS) cvodes_error("CVodeQuadSStolerances", flag);
        }
      }
    }

    // Adjoint sensitivity problem
    if (nrx1_>0) {
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

  int CvodesInterface::rhs(double t, N_Vector x, N_Vector xdot, void *user_data) {
    try {
      casadi_assert(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      const double* arg[] = {&t, NV_DATA_S(x), get_ptr(m->p)};
      double* res[] = {NV_DATA_S(xdot)};
      s.calc_function(m, "odeF", arg, res);
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
    if (nq1_>0) {
      N_VConst(0.0, m->q);
      flag = CVodeQuadReInit(m->mem, m->q);
      if (flag != CV_SUCCESS) cvodes_error("CVodeQuadReInit", flag);
    }

    // Re-initialize sensitivities
    if (ns_>0) {
      flag = CVodeSensReInit(m->mem, CV_SIMULTANEOUS, m->xzS);
      if (flag != CV_SUCCESS) cvodes_error("CVodeSensReInit", flag);
    }

    // Re-initialize backward integration
    if (nrx1_>0) {
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
      if (nrx1_>0) {
        // ... with taping
        int flag = CVodeF(m->mem, t, m->xz, &m->t, CV_NORMAL, &m->ncheck);
        if (flag!=CV_SUCCESS && flag!=CV_TSTOP_RETURN) cvodes_error("CVodeF", flag);
      } else {
        // ... without taping
        int flag = CVode(m->mem, t, m->xz, &m->t, CV_NORMAL);
        if (flag!=CV_SUCCESS && flag!=CV_TSTOP_RETURN) cvodes_error("CVode", flag);

        // Get sensitivities
        if (ns_>0) {
          double tret;
          int flag = CVodeGetSens(m->mem, &tret, m->xzS);
          if (flag!=CV_SUCCESS) cvodes_error("CVodeGetSens", flag);
        }
      }

      // Get quadratures
      if (nq1_>0) {
        double tret;
        int flag = CVodeGetQuad(m->mem, &tret, m->q);
        if (flag!=CV_SUCCESS) cvodes_error("CVodeGetQuad", flag);

        // sensitivities
        if (ns_>0) {
          int flag = CVodeGetQuadSens(m->mem, &tret, m->qS);
          if (flag!=CV_SUCCESS) cvodes_error("CVodeGetQuadSens", flag);
        }
      }
    }

    // Set function outputs
    casadi_copy(NV_DATA_S(m->xz), nx1_, x);
    for (int d=0; d<ns_; ++d) casadi_copy(NV_DATA_S(m->xzS[d]), nx1_, x+(d+1)*nx1_);
    casadi_copy(NV_DATA_S(m->q), nq1_, q);
    for (int d=0; d<ns_; ++d) casadi_copy(NV_DATA_S(m->qS[d]), nq1_, q+(d+1)*nq1_);

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
      if (nrq1_>0) {
        flag = CVodeGetQuadB(m->mem, m->whichB, &m->t, m->rq);
        if (flag!=CV_SUCCESS) cvodes_error("CVodeGetQuadB", flag);
      }
    }

    // Save outputs
    casadi_copy(NV_DATA_S(m->rxz), nrx1_, rx);
    casadi_copy(NV_DATA_S(m->rq), nrq1_, rq);

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

  void CvodesInterface::ehfun(int error_code, const char *module, const char *function,
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

  int CvodesInterface::rhsS(int Ns, double t, N_Vector x, N_Vector xdot, N_Vector *xS,
                            N_Vector *xdotS, void *user_data,
                            N_Vector tmp1, N_Vector tmp2) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;
      // Inputs and outputs
      const double** arg = m->arg;
      double** res = m->res;
      // All parameters
      const double* p = get_ptr(m->p);
      // Nondifferentiated inputs
      *arg++ = &t;
      *arg++ = NV_DATA_S(x);
      *arg++ = p;
      p += s.np1_;
      // Nondifferentiated outputs
      *arg++ = NV_DATA_S(xdot);
      // Forward seeds
      for (int d=0; d<s.ns_; ++d) {
        *arg++ = 0;
        *arg++ = NV_DATA_S(xS[d]);
        *arg++ = p;
        p += s.np1_;
      }
      // Forward sensitivities
      for (int d=0; d<s.ns_; ++d) {
        *res++ = NV_DATA_S(xdotS[d]);
      }
      // Evaluate
      s.calc_function(m, "odeFS");
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "fs failed: " << e.what() << endl;
      return 1;
    }
  }

  int CvodesInterface::rhsS1(int Ns, double t, N_Vector x, N_Vector xdot, int iS,
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

  int CvodesInterface::rhsQ(double t, N_Vector x, N_Vector qdot, void *user_data) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;
      const double* arg[] = {&t, NV_DATA_S(x), get_ptr(m->p)};
      double* res[] = {NV_DATA_S(qdot)};
      s.calc_function(m, "quadF", arg, res);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhsQ failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::rhsQS(int Ns, double t, N_Vector x, N_Vector *xS, N_Vector qdot,
                             N_Vector *qdotS, void *user_data,
                             N_Vector tmp1, N_Vector tmp2) {
    try {
      if (!user_data) {
        // SUNDIALS BUG!!!
        for (int i=0; i<Ns; ++i) N_VConst(0.0, qdotS[i]);
        return 0;
      }
      auto m = to_mem(user_data);
      auto& s = m->self;
      // Inputs and outputs
      const double** arg = m->arg;
      double** res = m->res;
      // All parameters
      const double* p = get_ptr(m->p);
      // Nondifferentiated inputs
      *arg++ = &t;
      *arg++ = NV_DATA_S(x);
      *arg++ = p;
      p += s.np1_;
      // Nondifferentiated outputs
      *arg++ = NV_DATA_S(qdot);
      // Forward seeds
      for (int d=0; d<s.ns_; ++d) {
        *arg++ = 0;
        *arg++ = NV_DATA_S(xS[d]);
        *arg++ = p;
        p += s.np1_;
      }
      // Forward sensitivities
      for (int d=0; d<s.ns_; ++d) {
        *res++ = NV_DATA_S(qdotS[d]);
      }
      // Evaluate
      s.calc_function(m, "quadFS");
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhsQS failed: " << e.what() << endl;;
      return 1;
    }
  }


  int CvodesInterface::rhsB(double t, N_Vector x, N_Vector rx, N_Vector rxdot,
                            void *user_data) {
    try {
      casadi_assert(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      const double* arg[] = {&t, NV_DATA_S(rx), get_ptr(m->rp),
        NV_DATA_S(x), get_ptr(m->p)};
      double* res[] = {NV_DATA_S(rxdot)};
      s.calc_function(m, "odeB", arg, res);

      // Negate (note definition of g)
      casadi_scal(s.nrx1_, -1., NV_DATA_S(rxdot));

      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhsB failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::rhsBS(double t, N_Vector x, N_Vector *xF, N_Vector xB,
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

  int CvodesInterface::rhsQB(double t, N_Vector x, N_Vector rx,
                             N_Vector rqdot, void *user_data) {
    try {
      casadi_assert(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      const double* arg[] = {&t, NV_DATA_S(rx), get_ptr(m->rp),
        NV_DATA_S(x), get_ptr(m->p)};
      double* res[] = {NV_DATA_S(rqdot)};
      s.calc_function(m, "quadB", arg, res);

      // Negate (note definition of g)
      casadi_scal(s.nrq1_, -1., NV_DATA_S(rqdot));

      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "rhsQB failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::jtimes(N_Vector v, N_Vector Jv, double t, N_Vector x,
                              N_Vector xdot, void *user_data, N_Vector tmp) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;
      const double* arg[] = {&t, NV_DATA_S(x), get_ptr(m->p), NV_DATA_S(v)};
      double* res[] = {NV_DATA_S(Jv)};
      s.calc_function(m, "jtimesF", arg, res);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "jtimes failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::jtimesB(N_Vector v, N_Vector Jv, double t, N_Vector x,
                               N_Vector rx, N_Vector rxdot, void *user_data ,
                               N_Vector tmpB) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;
      const double* arg[] = {&t, NV_DATA_S(x), get_ptr(m->p),
        NV_DATA_S(rx), get_ptr(m->rp), NV_DATA_S(v)};
      double* res[] = {NV_DATA_S(Jv)};
      s.calc_function(m, "jtimesB", arg, res);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "jtimes failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::djac(long N, double t, N_Vector x, N_Vector xdot, DlsMat Jac,
                            void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    try {
      casadi_assert(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      double one=1;
      const double* arg[] = {&t, NV_DATA_S(x), get_ptr(m->p), &one, 0};
      double* res[] = {m->jac};
      s.calc_function(m, "jacF", arg, res);

      // Save to Jac
      const Sparsity& sp = s.get_function("jacF").sparsity_out(0);
      const int* colind = sp.colind();
      int ncol = sp.size2();
      const int* row = sp.row();
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

  int CvodesInterface::djacB(long NeqB, double t, N_Vector x, N_Vector rx, N_Vector rxdot,
                             DlsMat JacB, void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                             N_Vector tmp3B) {
    try {
      casadi_assert(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      double minus_one = -1;
      const double* arg[] = {&t, NV_DATA_S(rx), get_ptr(m->rp),
            NV_DATA_S(x), get_ptr(m->p), &minus_one, 0};
      double* res[] = {m->jacB};
      s.calc_function(m, "jacB", arg, res);

      // Save to JacB
      const Sparsity& sp = s.get_function("jacB").sparsity_out(0);
      const int* colind = sp.colind();
      int ncol = sp.size2();
      const int* row = sp.row();
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

  int CvodesInterface::bjac(long N, long mupper, long mlower, double t, N_Vector x,
                            N_Vector xdot, DlsMat Jac, void *user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    try {
      casadi_assert(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      double one=1;
      const double* arg[] = {&t, NV_DATA_S(x), get_ptr(m->p), &one, 0};
      double* res[] = {m->jac};
      s.calc_function(m, "jacF", arg, res);

      // Save to Jac
      const Sparsity& sp = s.get_function("jacF").sparsity_out(0);
      const int* colind = sp.colind();
      int ncol = sp.size2();
      const int* row = sp.row();
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

  int CvodesInterface::bjacB(long NeqB, long mupperB, long mlowerB, double t, N_Vector x,
                             N_Vector rx, N_Vector rxdot, DlsMat JacB, void *user_data,
                             N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
    try {
      casadi_assert(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      double minus_one = -1;
      const double* arg[] = {&t, NV_DATA_S(rx), get_ptr(m->rp),
            NV_DATA_S(x), get_ptr(m->p), &minus_one, 0};
      double* res[] = {m->jacB};
      s.calc_function(m, "jacB", arg, res);

      // Save to JacB
      const Sparsity& sp = s.get_function("jacB").sparsity_out(0);
      const int* colind = sp.colind();
      int ncol = sp.size2();
      const int* row = sp.row();
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

  int CvodesInterface::psolve(double t, N_Vector x, N_Vector xdot, N_Vector r,
                              N_Vector z, double gamma, double delta, int lr,
                              void *user_data, N_Vector tmp) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;
      // Copy input to output, if necessary
      if (r!=z) {
        N_VScale(1.0, r, z);
      }

      // Solve the (possibly factorized) system
      const Function& linsol = s.get_function("linsolF");
      casadi_assert(linsol.nnz_out(0) == NV_LENGTH_S(z));
      linsol.linsol_solve(NV_DATA_S(z));

      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "psolve failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::psolveB(double t, N_Vector x, N_Vector xB, N_Vector xdotB,
                               N_Vector rvecB, N_Vector zvecB, double gammaB,
                               double deltaB, int lr, void *user_data, N_Vector tmpB) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;
      // Copy input to output, if necessary
      if (rvecB!=zvecB) {
        N_VScale(1.0, rvecB, zvecB);
      }

      // Solve the (possibly factorized) system
      const Function& linsolB = s.get_function("linsolB");
      casadi_assert(linsolB.nnz_out(0) == NV_LENGTH_S(zvecB));
      linsolB.linsol_solve(NV_DATA_S(zvecB), 1);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "psolveB failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::psetup(double t, N_Vector x, N_Vector xdot, booleantype jok,
                              booleantype *jcurPtr, double gamma, void *user_data,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;
      // Calculate Jacobian
      double d1 = -gamma, d2 = 1.;
      const double* arg[] = {&t, NV_DATA_S(x), get_ptr(m->p), &d1, &d2};
      double* res[] = {m->jac};
      s.calc_function(m, "jacF", arg, res);

      // Prepare the solution of the linear system (e.g. factorize)
      const Function& linsol = s.get_function("linsolF");
      linsol.setup(m->arg+LINSOL_NUM_IN, m->res+LINSOL_NUM_OUT, m->iw, m->w);
      linsol.linsol_factorize(m->jac);

      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "psetup failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::psetupB(double t, N_Vector x, N_Vector rx, N_Vector rxdot,
                               booleantype jokB, booleantype *jcurPtrB, double gammaB,
                               void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                               N_Vector tmp3B) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;
      // Calculate Jacobian
      double one=1;
      const double* arg[] = {&t, NV_DATA_S(rx), get_ptr(m->rp),
            NV_DATA_S(x), get_ptr(m->p), &gammaB, &one};
      double* res[] = {m->jacB};
      s.calc_function(m, "jacB", arg, res);

      // Prepare the solution of the linear system (e.g. factorize)
      const Function& linsolB = s.get_function("linsolB");
      linsolB.setup(m->arg+LINSOL_NUM_IN, m->res+LINSOL_NUM_OUT, m->iw, m->w);
      linsolB.linsol_factorize(m->jacB);

      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "psetupB failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::lsetup(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
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
      if (psetup(t, x, xdot, FALSE, jcurPtr,
                 gamma, static_cast<void*>(m), vtemp1, vtemp2, vtemp3)) return 1;

      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "lsetup failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::lsetupB(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
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
      if (psetupB(t, ca_mem->ca_ytmp, x, xdot, FALSE, jcurPtr,
                  gamma, static_cast<void*>(m), vtemp1, vtemp2, vtemp3)) return 1;

      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "lsetupB failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::lsolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                              N_Vector x, N_Vector xdot) {
    try {
      auto m = to_mem(cv_mem->cv_lmem);
      auto& s = m->self;

      // Current time
      double t = cv_mem->cv_tn;

      // Scaling factor before J
      double gamma = cv_mem->cv_gamma;

      // Accuracy
      double delta = 0.0;

      // Left/right preconditioner
      int lr = 1;

      // Call the preconditioner solve function (which solves the linear system)
      if (psolve(t, x, xdot, b, b, gamma, delta,
                 lr, static_cast<void*>(m), 0)) return 1;

      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "lsolve failed: " << e.what() << endl;;
      return 1;
    }
  }

  int CvodesInterface::lsolveB(CVodeMem cv_mem, N_Vector b, N_Vector weight,
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



      // Accuracy
      double delta = 0.0;

      // Left/right preconditioner
      int lr = 1;

      // Call the preconditioner solve function (which solves the linear system)
      if (psolveB(t, ca_mem->ca_ytmp, x, xdot, b, b, gamma, delta, lr,
                  static_cast<void*>(m), 0)) return 1;

      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "lsolveB failed: " << e.what() << endl;;
      return 1;
    }
  }

  void CvodesInterface::initDenseLinsol(CvodesMemory* m) const {
    int flag = CVDense(m->mem, nx1_);
    if (flag!=CV_SUCCESS) cvodes_error("CVDense", flag);
    if (exact_jacobian_) {
      flag = CVDlsSetDenseJacFn(m->mem, djac);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsSetDenseJacFn", flag);
    }
  }

  void CvodesInterface::initBandedLinsol(CvodesMemory* m) const {
    pair<int, int> bw = getBandwidth();
    int flag = CVBand(m->mem, nx1_, bw.first, bw.second);
    if (flag!=CV_SUCCESS) cvodes_error("CVBand", flag);
    if (exact_jacobian_) {
      flag = CVDlsSetBandJacFn(m->mem, bjac);
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
      flag = CVSpilsSetJacTimesVecFn(m->mem, jtimes);
      if (flag!=CV_SUCCESS) cvodes_error("CVSpilsSetJacTimesVecFn", flag);
    }

    // Add a preconditioner
    if (use_preconditioner_) {
      casadi_assert_message(has_function("jacF"), "No Jacobian function");
      casadi_assert_message(has_function("linsolF"), "No linear solver");
      flag = CVSpilsSetPreconditioner(m->mem, psetup, psolve);
      if (flag != CV_SUCCESS) cvodes_error("CVSpilsSetPreconditioner", flag);
    }
  }

  void CvodesInterface::initUserDefinedLinsol(CvodesMemory* m) const {
    casadi_assert_message(has_function("jacF"), "No Jacobian function");
    casadi_assert_message(has_function("linsolF"), "No linear solver");
    CVodeMem cv_mem = static_cast<CVodeMem>(m->mem);
    cv_mem->cv_lmem   = m;
    cv_mem->cv_lsetup = lsetup;
    cv_mem->cv_lsolve = lsolve;
    cv_mem->cv_setupNonNull = TRUE;
  }

  void CvodesInterface::initDenseLinsolB(CvodesMemory* m) const {
    int flag = CVDenseB(m->mem, m->whichB, nrx1_);
    if (flag!=CV_SUCCESS) cvodes_error("CVDenseB", flag);
    if (exact_jacobianB_) {
      flag = CVDlsSetDenseJacFnB(m->mem, m->whichB, djacB);
      if (flag!=CV_SUCCESS) cvodes_error("CVDlsSetDenseJacFnB", flag);
    }
  }

  void CvodesInterface::initBandedLinsolB(CvodesMemory* m) const {
    pair<int, int> bw = getBandwidthB();
    int flag = CVBandB(m->mem, m->whichB, nrx1_, bw.first, bw.second);
    if (flag!=CV_SUCCESS) cvodes_error("CVBandB", flag);

    if (exact_jacobianB_) {
      flag = CVDlsSetBandJacFnB(m->mem, m->whichB, bjacB);
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
      flag = CVSpilsSetJacTimesVecFnB(m->mem, m->whichB, jtimesB);
      if (flag!=CV_SUCCESS) cvodes_error("CVSpilsSetJacTimesVecFnB", flag);
    }

    // Add a preconditioner
    if (use_preconditionerB_) {
      casadi_assert_message(has_function("jacB"), "No Jacobian function");
      casadi_assert_message(has_function("linsolB"), "No linear solver");
      flag = CVSpilsSetPreconditionerB(m->mem, m->whichB, psetupB, psolveB);
      if (flag != CV_SUCCESS) cvodes_error("CVSpilsSetPreconditionerB", flag);
    }

  }

  void CvodesInterface::initUserDefinedLinsolB(CvodesMemory* m) const {
    casadi_assert_message(has_function("jacB"), "No Jacobian function");
    casadi_assert_message(has_function("linsolB"), "No linear solver");
    CVodeMem cv_mem = static_cast<CVodeMem>(m->mem);
    CVadjMem ca_mem = cv_mem->cv_adj_mem;
    CVodeBMem cvB_mem = ca_mem->cvB_mem;
    cvB_mem->cv_lmem   = m;
    cvB_mem->cv_mem->cv_lmem = m;
    cvB_mem->cv_mem->cv_lsetup = lsetupB;
    cvB_mem->cv_mem->cv_lsolve = lsolveB;
    cvB_mem->cv_mem->cv_setupNonNull = TRUE;
  }

  template<typename MatType>
  Function CvodesInterface::getJacF() {
    vector<MatType> a = MatType::get_input(oracle_);
    vector<MatType> r = oracle_(a);

    // Get the Jacobian in the Newton iteration
    MatType c_x = MatType::sym("c_x");
    MatType c_xdot = MatType::sym("c_xdot");
    MatType jac = c_x*MatType::jacobian(r[DE_ODE], a[DE_X]) + c_xdot*MatType::eye(nx1_);

    // Return generated function
    return Function("jacF", {a[DE_T], a[DE_X], a[DE_P], c_x, c_xdot}, {jac});
  }

  template<typename MatType>
  Function CvodesInterface::getJacB() {
    vector<MatType> a = MatType::get_input(oracle_);
    vector<MatType> r = oracle_(a);

    // Get the Jacobian in the Newton iteration
    MatType c_x = MatType::sym("c_x");
    MatType c_xdot = MatType::sym("c_xdot");
    MatType jac = c_x*MatType::jacobian(r[DE_RODE], a[DE_RX]) + c_xdot*MatType::eye(nrx1_);

    // return generated function
    return Function("jacB",
                    {a[DE_T], a[DE_RX], a[DE_RP], a[DE_X], a[DE_P], c_x, c_xdot}, {jac});
  }

  CvodesMemory::CvodesMemory(const CvodesInterface& s) : SundialsMemory(s), self(s) {
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
