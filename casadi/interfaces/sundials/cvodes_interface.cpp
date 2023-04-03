/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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
#include "casadi/core/casadi_misc.hpp"

#define THROWING(fcn, ...) \
cvodes_error(CASADI_STR(fcn), fcn(__VA_ARGS__))

namespace casadi {

extern "C"
int CASADI_INTEGRATOR_CVODES_EXPORT
casadi_register_integrator_cvodes(Integrator::Plugin* plugin) {
  plugin->creator = CvodesInterface::creator;
  plugin->name = "cvodes";
  plugin->doc = CvodesInterface::meta_doc.c_str();;
  plugin->version = CASADI_VERSION;
  plugin->options = &CvodesInterface::options_;
  plugin->deserialize = &CvodesInterface::deserialize;
  return 0;
}

extern "C"
void CASADI_INTEGRATOR_CVODES_EXPORT casadi_load_integrator_cvodes() {
  Integrator::registerPlugin(casadi_register_integrator_cvodes);
}

CvodesInterface::CvodesInterface(const std::string& name, const Function& dae,
    double t0, const std::vector<double>& tout) : SundialsInterface(name, dae, t0, tout) {
}

CvodesInterface::~CvodesInterface() {
  clear_mem();
}

const Options CvodesInterface::options_
= {{&SundialsInterface::options_},
    {{"linear_multistep_method",
      {OT_STRING,
      "Integrator scheme: BDF|adams"}},
    {"nonlinear_solver_iteration",
      {OT_STRING,
      "Nonlinear solver type: NEWTON|functional"}},
    {"min_step_size",
      {OT_DOUBLE,
      "Min step size [default: 0/0.0]"}},
    {"fsens_all_at_once",
      {OT_BOOL,
      "Calculate all right hand sides of the sensitivity equations at once"}}
    }
};

void CvodesInterface::init(const Dict& opts) {
  if (verbose_) casadi_message(name_ + "::init");

  // Initialize the base classes
  SundialsInterface::init(opts);

  // Default options
  std::string linear_multistep_method = "bdf";
  std::string nonlinear_solver_iteration = "newton";
  min_step_size_ = 0;

  // Read options
  for (auto&& op : opts) {
    if (op.first=="linear_multistep_method") {
      linear_multistep_method = op.second.to_string();
    } else if (op.first=="min_step_size") {
      min_step_size_ = op.second;
    } else if (op.first=="nonlinear_solver_iteration") {
      nonlinear_solver_iteration = op.second.to_string();
    }
  }

  // Algebraic variables not supported
  casadi_assert(nz_==0 && nrz_==0,
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

  // Misc
  alloc_w(nx_); // casadi_project
  alloc_w(nrx_); // casadi_project
}

int CvodesInterface::init_mem(void* mem) const {
  if (SundialsInterface::init_mem(mem)) return 1;
  auto m = to_mem(mem);

  // Create CVodes memory block
  m->mem = CVodeCreate(lmm_, iter_);
  casadi_assert(m->mem!=nullptr, "CVodeCreate: Creation failed");

  // Set error handler function
  THROWING(CVodeSetErrHandlerFn, m->mem, ehfun, m);

  // Set user data
  THROWING(CVodeSetUserData, m->mem, m);

  // Initialize CVodes
  double t0 = 0;
  THROWING(CVodeInit, m->mem, rhsF, t0, m->xz);

  // Set tolerances
  THROWING(CVodeSStolerances, m->mem, reltol_, abstol_);

  // Maximum number of steps
  THROWING(CVodeSetMaxNumSteps, m->mem, max_num_steps_);

  // Initial step size
  if (step0_!=0) THROWING(CVodeSetInitStep, m->mem, step0_);

  // Min step size
  if (min_step_size_!=0) THROWING(CVodeSetMinStep, m->mem, min_step_size_);

  // Max step size
  if (max_step_size_!=0) THROWING(CVodeSetMaxStep, m->mem, max_step_size_);

  // Maximum order of method
  if (max_order_) THROWING(CVodeSetMaxOrd, m->mem, max_order_);

  // Coeff. in the nonlinear convergence test
  if (nonlin_conv_coeff_!=0) THROWING(CVodeSetNonlinConvCoef, m->mem, nonlin_conv_coeff_);

  // attach a linear solver
  if (newton_scheme_==SD_DIRECT) {
    // Direct scheme
    CVodeMem cv_mem = static_cast<CVodeMem>(m->mem);
    cv_mem->cv_lmem   = m;
    cv_mem->cv_lsetup = lsetupF;
    cv_mem->cv_lsolve = lsolveF;
    cv_mem->cv_setupNonNull = TRUE;
  } else {
    // Iterative scheme
    casadi_int pretype = use_precon_ ? PREC_LEFT : PREC_NONE;
    switch (newton_scheme_) {
    case SD_DIRECT: casadi_assert_dev(0);
    case SD_GMRES: THROWING(CVSpgmr, m->mem, pretype, max_krylov_); break;
    case SD_BCGSTAB: THROWING(CVSpbcg, m->mem, pretype, max_krylov_); break;
    case SD_TFQMR: THROWING(CVSptfqmr, m->mem, pretype, max_krylov_); break;
    }
    THROWING(CVSpilsSetJacTimesVecFn, m->mem, jtimesF);
    if (use_precon_) THROWING(CVSpilsSetPreconditioner, m->mem, psetupF, psolveF);
  }

  // Quadrature equations
  if (nq_>0) {
    // Initialize quadratures in CVodes
    THROWING(CVodeQuadInit, m->mem, rhsQF, m->q);

    // Should the quadrature errors be used for step size control?
    if (quad_err_con_) {
      THROWING(CVodeSetQuadErrCon, m->mem, true);

      // Quadrature error tolerances
      // TODO(Joel): vector absolute tolerances
      THROWING(CVodeQuadSStolerances, m->mem, reltol_, abstol_);
    }
  }

  // Initialize adjoint sensitivities
  if (nrx_>0) {
    casadi_int interpType = interp_ == SD_HERMITE ? CV_HERMITE : CV_POLYNOMIAL;
    THROWING(CVodeAdjInit, m->mem, steps_per_checkpoint_, interpType);
  }

  m->first_callB = true;
  return 0;
}

int CvodesInterface::rhsF(double t, N_Vector x, N_Vector xdot, void *user_data) {
  try {
    casadi_assert_dev(user_data);
    auto m = to_mem(user_data);
    auto& s = m->self;
    if (s.calc_daeF(m, t, NV_DATA_S(x), nullptr, NV_DATA_S(xdot), nullptr)) return 1;
    return 0;
  } catch(std::exception& e) { // non-recoverable error
    uerr() << "rhs failed: " << e.what() << std::endl;
    return -1;
  }
}

void CvodesInterface::reset(IntegratorMemory* mem,
    const double* x, const double* z, const double* _p) const {
  if (verbose_) casadi_message(name_ + "::reset");
  auto m = to_mem(mem);

  // Reset the base classes
  SundialsInterface::reset(mem, x, z, _p);

  // Re-initialize
  THROWING(CVodeReInit, m->mem, m->t, m->xz);

  // Re-initialize quadratures
  if (nq_ > 0) {
    N_VConst(0.0, m->q);
    THROWING(CVodeQuadReInit, m->mem, m->q);
  }

  // Re-initialize backward integration
  if (nrx_ > 0) {
    THROWING(CVodeAdjReInit, m->mem);
  }
}

void CvodesInterface::advance(IntegratorMemory* mem,
    const double* u, double* x, double* z, double* q) const {
  auto m = to_mem(mem);

  // Set controls
  casadi_copy(u, nu_, m->u);

  // Do not integrate past change in input signals or past the end
  THROWING(CVodeSetStopTime, m->mem, m->t_stop);

  // Integrate, unless already at desired time
  const double ttol = 1e-9;
  if (fabs(m->t - m->t_next) >= ttol) {
    // Integrate forward ...
    double tret = m->t;
    if (nrx_>0) {
      // ... with taping
      THROWING(CVodeF, m->mem, m->t_next, m->xz, &tret, CV_NORMAL, &m->ncheck);
    } else {
      // ... without taping
      THROWING(CVode, m->mem, m->t_next, m->xz, &tret, CV_NORMAL);
    }

    // Get quadratures
    if (nq_ > 0) {
      THROWING(CVodeGetQuad, m->mem, &tret, m->q);
    }
  }

  // Set function outputs
  casadi_copy(NV_DATA_S(m->xz), nx_, x);
  casadi_copy(NV_DATA_S(m->q), nq_, q);

  // Get stats
  THROWING(CVodeGetIntegratorStats, m->mem, &m->nsteps, &m->nfevals, &m->nlinsetups,
            &m->netfails, &m->qlast, &m->qcur, &m->hinused,
            &m->hlast, &m->hcur, &m->tcur);
  THROWING(CVodeGetNonlinSolvStats, m->mem, &m->nniters, &m->nncfails);
}

void CvodesInterface::resetB(IntegratorMemory* mem,
    const double* rx, const double* rz, const double* rp) const {
  auto m = to_mem(mem);

  // Reset the base classes
  SundialsInterface::resetB(mem, rx, rz, rp);

  if (m->first_callB) {
    // Create backward problem
    THROWING(CVodeCreateB, m->mem, lmm_, iter_, &m->whichB);
    THROWING(CVodeInitB, m->mem, m->whichB, rhsB, m->t, m->rxz);
    THROWING(CVodeSStolerancesB, m->mem, m->whichB, reltol_, abstol_);
    THROWING(CVodeSetUserDataB, m->mem, m->whichB, m);
    if (newton_scheme_==SD_DIRECT) {
      // Direct scheme
      CVodeMem cv_mem = static_cast<CVodeMem>(m->mem);
      CVadjMem ca_mem = cv_mem->cv_adj_mem;
      CVodeBMem cvB_mem = ca_mem->cvB_mem;
      cvB_mem->cv_lmem = m;
      cvB_mem->cv_mem->cv_lmem = m;
      cvB_mem->cv_mem->cv_lsetup = lsetupB;
      cvB_mem->cv_mem->cv_lsolve = lsolveB;
      cvB_mem->cv_mem->cv_setupNonNull = TRUE;
    } else {
      // Iterative scheme
      casadi_int pretype = use_precon_ ? PREC_LEFT : PREC_NONE;
      switch (newton_scheme_) {
      case SD_DIRECT: casadi_assert_dev(0);
      case SD_GMRES: THROWING(CVSpgmrB, m->mem, m->whichB, pretype, max_krylov_); break;
      case SD_BCGSTAB: THROWING(CVSpbcgB, m->mem, m->whichB, pretype, max_krylov_); break;
      case SD_TFQMR: THROWING(CVSptfqmrB, m->mem, m->whichB, pretype, max_krylov_); break;
      }
      THROWING(CVSpilsSetJacTimesVecFnB, m->mem, m->whichB, jtimesB);
      if (use_precon_) THROWING(CVSpilsSetPreconditionerB, m->mem, m->whichB, psetupB, psolveB);
    }

    // Quadratures for the backward problem
    THROWING(CVodeQuadInitB, m->mem, m->whichB, rhsQB, m->ruq);
    if (quad_err_con_) {
      THROWING(CVodeSetQuadErrConB, m->mem, m->whichB, true);
      THROWING(CVodeQuadSStolerancesB, m->mem, m->whichB, reltol_, abstol_);
    }

    // Mark initialized
    m->first_callB = false;
  } else {
    THROWING(CVodeReInitB, m->mem, m->whichB, m->t, m->rxz);
    THROWING(CVodeQuadReInitB, m->mem, m->whichB, m->ruq);
  }
}

void CvodesInterface::impulseB(IntegratorMemory* mem,
    const double* rx, const double* rz, const double* rp) const {
  auto m = to_mem(mem);

  // Call method in base class
  SundialsInterface::impulseB(mem, rx, rz, rp);

  // Reinitialize solver
  THROWING(CVodeReInitB, m->mem, m->whichB, m->t, m->rxz);
  THROWING(CVodeQuadReInitB, m->mem, m->whichB, m->ruq);
}

void CvodesInterface::retreat(IntegratorMemory* mem, const double* u,
    double* rx, double* rz, double* rq, double* uq) const {
  auto m = to_mem(mem);

  // Set controls
  casadi_copy(u, nu_, m->u);

  // Integrate, unless already at desired time
  if (m->t_next < m->t) {
    THROWING(CVodeB, m->mem, m->t_next, CV_NORMAL);
    double tret;
    THROWING(CVodeGetB, m->mem, m->whichB, &tret, m->rxz);
    if (nrq_ > 0 || nuq_ > 0) {
      THROWING(CVodeGetQuadB, m->mem, m->whichB, &tret, m->ruq);
    }
  }

  // Save outputs
  casadi_copy(NV_DATA_S(m->rxz), nrx_, rx);
  casadi_copy(NV_DATA_S(m->ruq), nrq_, rq);
  casadi_copy(NV_DATA_S(m->ruq) + nrq_, nuq_, uq);

  // Get stats
  CVodeMem cv_mem = static_cast<CVodeMem>(m->mem);
  CVadjMem ca_mem = cv_mem->cv_adj_mem;
  CVodeBMem cvB_mem = ca_mem->cvB_mem;
  THROWING(CVodeGetIntegratorStats, cvB_mem->cv_mem, &m->nstepsB,
          &m->nfevalsB, &m->nlinsetupsB, &m->netfailsB, &m->qlastB,
          &m->qcurB, &m->hinusedB, &m->hlastB, &m->hcurB, &m->tcurB);
  THROWING(CVodeGetNonlinSolvStats, cvB_mem->cv_mem, &m->nnitersB, &m->nncfailsB);
}

void CvodesInterface::cvodes_error(const char* module, int flag) {
  // Successfull return or warning
  if (flag>=CV_SUCCESS) return;
  // Construct error message
  char* flagname = CVodeGetReturnFlagName(flag);
  std::stringstream ss;
  ss << module << " returned \"" << flagname << "\". Consult CVODES documentation.";
  free(flagname);  // NOLINT
  casadi_error(ss.str());
}

void CvodesInterface::ehfun(int error_code, const char *module, const char *function,
    char *msg, void *user_data) {
  try {
    casadi_assert_dev(user_data);
    auto m = to_mem(user_data);
    auto& s = m->self;
    if (!s.disable_internal_warnings_) {
      uerr() << msg << std::endl;
    }
  } catch(std::exception& e) {
    uerr() << "ehfun failed: " << e.what() << std::endl;
  }
}

int CvodesInterface::rhsQF(double t, N_Vector x, N_Vector qdot, void *user_data) {
  try {
    auto m = to_mem(user_data);
    auto& s = m->self;
    if (s.calc_quadF(m, t, NV_DATA_S(x), nullptr, NV_DATA_S(qdot))) return 1;

    return 0;
  } catch(std::exception& e) { // non-recoverable error
    uerr() << "rhsQ failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesInterface::rhsB(double t, N_Vector x, N_Vector rx, N_Vector rxdot, void *user_data) {
  try {
    casadi_assert_dev(user_data);
    auto m = to_mem(user_data);
    auto& s = m->self;
    if (s.calc_daeB(m, t, NV_DATA_S(x), nullptr, NV_DATA_S(rx), nullptr,
      NV_DATA_S(rxdot), nullptr)) return 1;
    // Negate (note definition of g)
    casadi_scal(s.nrx_, -1., NV_DATA_S(rxdot));
    return 0;
  } catch(std::exception& e) { // non-recoverable error
    uerr() << "rhsB failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesInterface::rhsQB(double t, N_Vector x, N_Vector rx, N_Vector ruqdot, void *user_data) {
  try {
    casadi_assert_dev(user_data);
    auto m = to_mem(user_data);
    auto& s = m->self;
    if (s.calc_quadB(m, t, NV_DATA_S(x), nullptr, NV_DATA_S(rx), nullptr,
      NV_DATA_S(ruqdot), NV_DATA_S(ruqdot) + s.nrq_)) return 1;

    // Negate (note definition of g)
    casadi_scal((s.nrq_ + s.nuq_), -1., NV_DATA_S(ruqdot));

    return 0;
  } catch(std::exception& e) { // non-recoverable error
    uerr() << "rhsQB failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesInterface::jtimesF(N_Vector v, N_Vector Jv, double t, N_Vector x,
    N_Vector xdot, void *user_data, N_Vector tmp) {
  try {
    auto m = to_mem(user_data);
    auto& s = m->self;
    if (s.calc_jtimesF(m, t, NV_DATA_S(x), nullptr, NV_DATA_S(v), nullptr,
      NV_DATA_S(Jv), nullptr)) return 1;
    return 0;
  } catch(std::exception& e) { // non-recoverable error
    uerr() << "jtimesF failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesInterface::jtimesB(N_Vector v, N_Vector Jv, double t, N_Vector x,
    N_Vector rx, N_Vector rxdot, void *user_data, N_Vector tmpB) {
  try {
    auto m = to_mem(user_data);
    auto& s = m->self;
    if (s.calc_jtimesB(m, t, NV_DATA_S(x), nullptr, NV_DATA_S(rx), nullptr,
      NV_DATA_S(v), nullptr, NV_DATA_S(Jv), nullptr)) return 1;
    return 0;
  } catch(std::exception& e) { // non-recoverable error
    uerr() << "jtimesB failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesInterface::psolveF(double t, N_Vector x, N_Vector xdot, N_Vector r,
    N_Vector z, double gamma, double delta, int lr, void *user_data, N_Vector tmp) {
  try {
    auto m = to_mem(user_data);
    auto& s = m->self;

    // Get right-hand sides in m->v1
    double* v = NV_DATA_S(r);
    casadi_copy(v, s.nx_, m->v1);

    // Solve for undifferentiated right-hand-side, save to output
    if (s.linsolF_.solve(m->jacF, m->v1, 1, false, m->mem_linsolF)) return 1;
    v = NV_DATA_S(z); // possibly different from r
    casadi_copy(m->v1, s.nx1_, v);

    // Sensitivity equations
    if (s.nfwd_ > 0) {
      // Second order correction
      if (s.second_order_correction_) {
        // The outputs will double as seeds for jtimesF
        casadi_clear(v + s.nx1_, s.nx_ - s.nx1_);
        if (s.calc_jtimesF(m, t, NV_DATA_S(x), nullptr, v, nullptr, m->v2, nullptr)) return 1;

        // Subtract m->v2 from m->v1, scaled with -gamma
        casadi_axpy(s.nx_ - s.nx1_, m->gamma, m->v2 + s.nx1_, m->v1 + s.nx1_);
      }

      // Solve for sensitivity right-hand-sides
      if (s.linsolF_.solve(m->jacF, m->v1 + s.nx1_, s.nfwd_, false, m->mem_linsolF)) return 1;

      // Save to output, reordered
      casadi_copy(m->v1 + s.nx1_, s.nx_ - s.nx1_, v + s.nx1_);
    }

    return 0;
  } catch(std::exception& e) { // non-recoverable error
    uerr() << "psolve failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesInterface::psolveB(double t, N_Vector x, N_Vector xB, N_Vector xdotB, N_Vector rvecB,
    N_Vector zvecB, double gammaB, double deltaB, int lr, void *user_data, N_Vector tmpB) {
  try {
    auto m = to_mem(user_data);
    auto& s = m->self;

    // Get right-hand sides in m->v1
    double* v = NV_DATA_S(rvecB);
    casadi_copy(v, s.nrx_, m->v1);

    // Solve for undifferentiated right-hand-side, save to output
    if (s.linsolF_.solve(m->jacF, m->v1, s.nadj_, true, m->mem_linsolF)) return 1;
    v = NV_DATA_S(zvecB); // possibly different from rvecB
    casadi_copy(m->v1, s.nrx1_ * s.nadj_, v);

    // Sensitivity equations
    if (s.nfwd_ > 0) {
      // Second order correction
      if (s.second_order_correction_) {
        // The outputs will double as seeds for jtimesB
        casadi_clear(v + s.nrx1_ * s.nadj_, s.nrx_ - s.nrx1_ * s.nadj_);
        if (s.calc_jtimesB(m, t, NV_DATA_S(x), nullptr, NV_DATA_S(xB),
          nullptr, v, nullptr, m->v2, nullptr)) return 1;

        // Subtract m->v2 from m->v1, scaled with gammaB
        casadi_axpy(s.nrx_ - s.nrx1_ * s.nadj_, -m->gammaB, m->v2 + s.nrx1_ * s.nadj_,
          m->v1 + s.nrx1_ * s.nadj_);
      }

      // Solve for sensitivity right-hand-sides
      if (s.linsolF_.solve(m->jacF, m->v1 + s.nx1_, s.nadj_ * s.nfwd_,
        true, m->mem_linsolF)) return 1;

      // Save to output, reordered
      casadi_copy(m->v1 + s.nx1_, s.nx_ - s.nx1_, v + s.nx1_);
    }

    return 0;
  } catch(std::exception& e) { // non-recoverable error
    uerr() << "psolveB failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesInterface::psetupF(double t, N_Vector x, N_Vector xdot, booleantype jok,
    booleantype *jcurPtr, double gamma, void *user_data,
    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  try {
    auto m = to_mem(user_data);
    auto& s = m->self;
    // Store gamma for later
    m->gamma = gamma;

    // Sparsity patterns
    const Sparsity& sp_jac_ode_x = s.get_function("jacF").sparsity_out(0);
    const Sparsity& sp_jacF = s.linsolF_.sparsity();

    // Offset for storing the sparser Jacobian, to allow overwriting entries
    casadi_int jac_offset = sp_jacF.nnz() - sp_jac_ode_x.nnz();

    // Calculate Jacobian
    if (s.calc_jacF(m, t, NV_DATA_S(x), nullptr,
      m->jacF + jac_offset, nullptr, nullptr, nullptr)) return 1;

    // Project to expected sparsity pattern (with diagonal)
    casadi_project(m->jacF + jac_offset, sp_jac_ode_x, m->jacF, sp_jacF, m->w);

    // Scale and shift diagonal
    const casadi_int *colind = sp_jacF.colind(), *row = sp_jacF.row();
    for (casadi_int c = 0; c < sp_jacF.size2(); ++c) {
      for (casadi_int k = colind[c]; k < colind[c + 1]; ++k) {
        casadi_int r = row[k];
        // Scale Jacobian
        m->jacF[k] *= -gamma;
        // Add contribution to diagonal
        if (r == c) m->jacF[k] += 1;
      }
    }

    // Jacobian is now current
    if (jcurPtr) *jcurPtr = 1;

    // Prepare the solution of the linear system (e.g. factorize)
    if (s.linsolF_.nfact(m->jacF, m->mem_linsolF)) return 1;

    return 0;
  } catch(std::exception& e) { // non-recoverable error
    uerr() << "psetup failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesInterface::psetupB(double t, N_Vector x, N_Vector rx, N_Vector rxdot,
    booleantype jokB, booleantype *jcurPtrB, double gammaB,
    void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
  try {
    auto m = to_mem(user_data);
    // Store gamma for later
    m->gammaB = gammaB;
    // We use the same linear solver for the forward problem as for the backward problem
    return psetupF(t, x, nullptr, jokB, jcurPtrB, -gammaB, user_data, tmp1B, tmp2B, tmp3B);

  } catch(std::exception& e) { // non-recoverable error
    uerr() << "psetupB failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesInterface::lsetupF(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
    booleantype *jcurPtr, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3) {
  try {
    auto m = to_mem(cv_mem->cv_lmem);

    // Call the preconditioner setup function (which sets up the linear solver)
    return psetupF(cv_mem->cv_tn, x, xdot, FALSE, jcurPtr,
      cv_mem->cv_gamma, static_cast<void*>(m), vtemp1, vtemp2, vtemp3);

  } catch(std::exception& e) { // non-recoverable error
    uerr() << "lsetup failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesInterface::lsetupB(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
    booleantype *jcurPtr, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3) {
  try {
    auto m = to_mem(cv_mem->cv_lmem);
    CVadjMem ca_mem;
    //CVodeBMem cvB_mem;

    // Current time
    double t = cv_mem->cv_tn; // TODO(Joel): is this correct?
    double gamma = cv_mem->cv_gamma;

    cv_mem = static_cast<CVodeMem>(cv_mem->cv_user_data);
    ca_mem = cv_mem->cv_adj_mem;
    //cvB_mem = ca_mem->ca_bckpbCrt;

    // Get FORWARD solution from interpolation.
    int flag = ca_mem->ca_IMget(cv_mem, t, ca_mem->ca_ytmp, nullptr);
    if (flag != CV_SUCCESS) casadi_error("Could not interpolate forward states");

    // Call the preconditioner setup function (which sets up the linear solver)
    return psetupB(t, ca_mem->ca_ytmp, x, xdot, FALSE, jcurPtr,
      gamma, static_cast<void*>(m), vtemp1, vtemp2, vtemp3);

  } catch(std::exception& e) { // non-recoverable error
    uerr() << "lsetupB failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesInterface::lsolveF(CVodeMem cv_mem, N_Vector b, N_Vector weight,
    N_Vector x, N_Vector xdot) {
  try {
    auto m = to_mem(cv_mem->cv_lmem);
    //auto& s = m->self;

    // Current time
    double t = cv_mem->cv_tn;

    // Scaling factor before J
    double gamma = cv_mem->cv_gamma;

    // Accuracy
    double delta = 0.0;

    // Left/right preconditioner
    casadi_int lr = 1;

    // Call the preconditioner solve function (which solves the linear system)
    return psolveF(t, x, xdot, b, b, gamma, delta, lr, static_cast<void*>(m), nullptr);

  } catch(std::exception& e) { // non-recoverable error
    uerr() << "lsolveF failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesInterface::lsolveB(CVodeMem cv_mem, N_Vector b, N_Vector weight,
    N_Vector x, N_Vector xdot) {
  try {
    auto m = to_mem(cv_mem->cv_lmem);
    CVadjMem ca_mem;
    //CVodeBMem cvB_mem;

    // Current time
    double t = cv_mem->cv_tn; // TODO(Joel): is this correct?
    double gamma = cv_mem->cv_gamma;

    cv_mem = static_cast<CVodeMem>(cv_mem->cv_user_data);

    ca_mem = cv_mem->cv_adj_mem;
    //cvB_mem = ca_mem->ca_bckpbCrt;

    // Get FORWARD solution from interpolation.
    int flag = ca_mem->ca_IMget(cv_mem, t, ca_mem->ca_ytmp, nullptr);
    if (flag != CV_SUCCESS) casadi_error("Could not interpolate forward states");

    // Accuracy
    double delta = 0.0;

    // Left/right preconditioner
    int lr = 1;

    // Call the preconditioner solve function (which solves the linear system)
    return psolveB(t, ca_mem->ca_ytmp, x, xdot, b, b, gamma, delta, lr,
      static_cast<void*>(m), nullptr);

  } catch(std::exception& e) { // non-recoverable error
    uerr() << "lsolveB failed: " << e.what() << std::endl;
    return -1;
  }
}

CvodesMemory::CvodesMemory(const CvodesInterface& s) : self(s) {
  this->mem = nullptr;

  // Reset checkpoints counter
  this->ncheck = 0;
}

CvodesMemory::~CvodesMemory() {
  if (this->mem_linsolF >= 0) self.linsolF_.release(this->mem_linsolF);
  if (this->mem) CVodeFree(&this->mem);
}

CvodesInterface::CvodesInterface(DeserializingStream& s) : SundialsInterface(s) {
  int version = s.version("CvodesInterface", 1, 2);
  s.unpack("CvodesInterface::lmm", lmm_);
  s.unpack("CvodesInterface::iter", iter_);

  if (version>=2) {
    s.unpack("CvodesInterface::min_step_size", min_step_size_);
  } else {
    min_step_size_ = 0;
  }
}

void CvodesInterface::serialize_body(SerializingStream &s) const {
  SundialsInterface::serialize_body(s);
  s.version("CvodesInterface", 2);

  s.pack("CvodesInterface::lmm", lmm_);
  s.pack("CvodesInterface::iter", iter_);
  s.pack("CvodesInterface::min_step_size", min_step_size_);
}

} // namespace casadi
