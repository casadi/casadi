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


#include "cvodes_simulator.hpp"
#include "casadi/core/casadi_misc.hpp"

#define THROWING(fcn, ...) \
cvodes_error(CASADI_STR(fcn), fcn(__VA_ARGS__))

namespace casadi {

extern "C"
int CASADI_SIMULATOR_CVODES_EXPORT
casadi_register_simulator_cvodes(Simulator::Plugin* plugin) {
  plugin->creator = CvodesSimulator::creator;
  plugin->name = "cvodes";
  plugin->doc = CvodesSimulator::meta_doc.c_str();;
  plugin->version = CASADI_VERSION;
  plugin->options = &CvodesSimulator::options_;
  return 0;
}

extern "C"
void CASADI_SIMULATOR_CVODES_EXPORT casadi_load_simulator_cvodes() {
  Simulator::registerPlugin(casadi_register_simulator_cvodes);
}

CvodesSimulator::CvodesSimulator(const std::string& name, const Function& dae,
  const std::vector<double>& grid)
  : SundialsSimulator(name, dae, grid) {
}

CvodesSimulator::~CvodesSimulator() {
  clear_mem();
}

const Options CvodesSimulator::options_
= {{&SundialsSimulator::options_},
   {{"linear_multistep_method",
     {OT_STRING,
      "Simulator scheme: BDF|adams"}},
    {"nonlinear_solver_iteration",
     {OT_STRING,
      "Nonlinear solver type: NEWTON|functional"}},
    {"min_step_size",
     {OT_DOUBLE,
      "Min step size [default: 0/0.0]"}}
   }
};

void CvodesSimulator::init(const Dict& opts) {
  if (verbose_) casadi_message(name_ + "::init");

  // Initialize the base classes
  SundialsSimulator::init(opts);

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
  casadi_assert(nz_==0, "CVODES does not support algebraic variables");

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
  if (newton_scheme_ != NewtonScheme::DIRECT) {
    set_function(oracle_.forward(1), "jtimes");
  }
}

int CvodesSimulator::init_mem(void* mem) const {
  if (SundialsSimulator::init_mem(mem)) return 1;
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
  THROWING(CVodeInit, m->mem, rhs, t0, m->xz);

  // Absolute tolerances for each component
  if (m->abstolv) {
    casadi_copy(get_ptr(nom_x_), nx_, NV_DATA_S(m->abstolv));
    casadi_copy(get_ptr(nom_z_), nz_, NV_DATA_S(m->abstolv) + nx_);
    casadi_scal(nx_+nz_, abstol_, NV_DATA_S(m->abstolv));
  }

  // Set tolerances
  if (scale_abstol_) {
    THROWING(CVodeSVtolerances, m->mem, reltol_, m->abstolv);
  } else {
    THROWING(CVodeSStolerances, m->mem, reltol_, abstol_);
  }

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
  if (newton_scheme_ == NewtonScheme::DIRECT) {
    // Direct scheme
    CVodeMem cv_mem = static_cast<CVodeMem>(m->mem);
    cv_mem->cv_lmem   = m;
    cv_mem->cv_lsetup = lsetup;
    cv_mem->cv_lsolve = lsolve;
    cv_mem->cv_setupNonNull = TRUE;
  } else {
    // Iterative scheme
    casadi_int pretype = use_precon_ ? PREC_LEFT : PREC_NONE;
    switch (newton_scheme_) {
      case NewtonScheme::DIRECT: casadi_assert_dev(0);
      case NewtonScheme::GMRES: THROWING(CVSpgmr, m->mem, pretype, max_krylov_); break;
      case NewtonScheme::BCGSTAB: THROWING(CVSpbcg, m->mem, pretype, max_krylov_); break;
      case NewtonScheme::TFQMR: THROWING(CVSptfqmr, m->mem, pretype, max_krylov_); break;
      default: casadi_error("No such Newton scheme");
    }
    THROWING(CVSpilsSetJacTimesVecFn, m->mem, jtimes);
    if (use_precon_) THROWING(CVSpilsSetPreconditioner, m->mem, psetup, psolve);
  }

  // Forward sensitivities
  if (nfwd_ > 0) {
    // Initialize sensitivity problem
    int ism = CV_SIMULTANEOUS;  // or CV_STAGGERED
    THROWING(CVodeSensInit, m->mem, nfwd_, ism, sens_rhs, get_ptr(m->fwd_xz));

    // Set tolerances
    if (scale_abstol_) {
      THROWING(CVodeSensSVtolerances, m->mem, reltol_, get_ptr(m->abstolvS));
    } else {
      casadi_assert(!abstolS_.empty(), "here");
      THROWING(CVodeSensSStolerances, m->mem, reltol_, const_cast<double*>(get_ptr(abstolS_)));
    }
  }

  return 0;
}

int CvodesSimulator::rhs(double t, N_Vector x, N_Vector xdot, void *user_data) {
  try {
    casadi_assert_dev(user_data);
    auto m = to_mem(user_data);
    auto& s = m->self;
    std::fill_n(m->arg, enum_traits<DynIn>::n_enum, nullptr);
    m->arg[DYN_T] = &t;
    m->arg[DYN_X] = NV_DATA_S(x);
    m->arg[DYN_U] = m->u;
    m->arg[DYN_P] = m->p;
    std::fill_n(m->res, enum_traits<DynOut>::n_enum, nullptr);
    m->res[DYN_ODE] = NV_DATA_S(xdot);
    s.calc_function(m, "dae");
    return 0;
  } catch(int flag) { // recoverable error
    return flag;
  } catch(std::exception& e) { // non-recoverable error
    uerr() << "rhs failed: " << e.what() << std::endl;
    return -1;
  }
}

void CvodesSimulator::reset(SimulatorMemory* mem) const {
  if (verbose_) casadi_message(name_ + "::reset");
  auto m = to_mem(mem);
  // Reset the base classes
  SundialsSimulator::reset(mem);
  // Re-initialize
  THROWING(CVodeReInit, m->mem, m->t, m->xz);
  // Forward sensitivities
  if (nfwd_ > 0) {
    int ism = CV_SIMULTANEOUS;  // or CV_STAGGERED
    THROWING(CVodeSensReInit, m->mem, ism, get_ptr(m->fwd_xz));
  }
}

void CvodesSimulator::advance(SimulatorMemory* mem, double t, double t_stop) const {
  auto m = to_mem(mem);
  // Do not integrate past change in input signals or past the end
  THROWING(CVodeSetStopTime, m->mem, t_stop);
  // Integrate
  THROWING(CVode, m->mem, t, m->xz, &m->t, CV_NORMAL);
  casadi_copy(NV_DATA_S(m->xz), nx_, m->xk);
  //Get forward sensitivities
  if (nfwd_ > 0) {
    THROWING(CVodeGetSens, m->mem, &m->t, get_ptr(m->fwd_xz));
    for (size_t i = 0; i < nfwd_; ++i) {
      casadi_copy(NV_DATA_S(m->fwd_xz[i]), nx_, m->fwd_xk + i * nx_);
    }
  }
  // Get stats
  THROWING(CVodeGetIntegratorStats, m->mem, &m->nsteps, &m->nfevals, &m->nlinsetups,
    &m->netfails, &m->qlast, &m->qcur, &m->hinused,
    &m->hlast, &m->hcur, &m->tcur);
  THROWING(CVodeGetNonlinSolvStats, m->mem, &m->nniters, &m->nncfails);
}

void CvodesSimulator::cvodes_error(const char* module, int flag) {
  // Successfull return or warning
  if (flag>=CV_SUCCESS) return;
  // Construct error message
  char* flagname = CVodeGetReturnFlagName(flag);
  std::stringstream ss;
  ss << module << " returned \"" << flagname << "\". Consult CVODES documentation.";
  free(flagname); // NOLINT
  casadi_error(ss.str());
}

void CvodesSimulator::ehfun(int error_code, const char *module, const char *function,
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

int CvodesSimulator::jtimes(N_Vector v, N_Vector Jv, double t, N_Vector x,
                            N_Vector xdot, void *user_data, N_Vector tmp) {
  try {
    auto m = to_mem(user_data);
    auto& s = m->self;
    // Set input and output buffers
    std::fill_n(m->arg, enum_traits<DynIn>::n_enum + enum_traits<DynOut>::n_enum
      + enum_traits<DynIn>::n_enum, nullptr);
    m->arg[DYN_T] = &t;  // t
    m->arg[DYN_X] = NV_DATA_S(x);  // x
    m->arg[DYN_U] = m->u;  // u
    m->arg[DYN_P] = m->p;  // p
    m->arg[enum_traits<DynIn>::n_enum + DYN_ODE] = NV_DATA_S(xdot);  // ode
    m->arg[enum_traits<DynIn>::n_enum
      + enum_traits<DynOut>::n_enum + DYN_X] = NV_DATA_S(v);  // fwd:x
    std::fill_n(m->res, enum_traits<DynOut>::n_enum, nullptr);
    m->res[DYN_ODE] = NV_DATA_S(Jv);  // fwd:ode
    // Evaluate
    s.calc_function(m, "jtimes");
    return 0;
  } catch(casadi_int flag) { // recoverable error
    return flag;
  } catch(std::exception& e) { // non-recoverable error
    uerr() << "jtimes failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesSimulator::psolve(double t, N_Vector x, N_Vector xdot, N_Vector r,
                            N_Vector z, double gamma, double delta, int lr,
                            void *user_data, N_Vector tmp) {
  try {
    auto m = to_mem(user_data);
    auto& s = m->self;

    // Get right-hand sides in m->v1
    double* v = NV_DATA_S(r);
    casadi_copy(v, s.nx_, m->v1);

    // Solve for undifferentiated right-hand-side, save to output
    if (s.linsolF_.solve(m->lin_nz, m->v1, 1, false, m->mem_linsolF))
      casadi_error("Linear system solve failed");
    v = NV_DATA_S(z); // possibly different from r
    casadi_copy(m->v1, s.nx_, v);

    return 0;
  } catch(int flag) { // recoverable error
    return flag;
  } catch(std::exception& e) { // non-recoverable error
    uerr() << "psolve failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesSimulator::psetup(double t, N_Vector x, N_Vector xdot, booleantype jok,
                            booleantype *jcurPtr, double gamma, void *user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  try {
    auto m = to_mem(user_data);
    auto& s = m->self;
    // Store gamma for later
    m->gamma = gamma;

    // Calculate Jacobian
    if (s.calculate_jac(m, t, x, xdot)) casadi_error("'jac' calculation failed");

    // Jacobian is now current
    *jcurPtr = 1;

    // Copy Jacobian and scale with -gamma
    casadi_copy(m->jac_x, s.jac_x_sp_.nnz(), m->lin_nz);
    casadi_scal(s.jac_x_sp_.nnz(), -gamma, m->lin_nz);

    // Add diagonal contribution, project to correct sparsity
    add_diag(s.jac_x_sp_, m->lin_nz, 1., s.lin_sp_, NV_DATA_S(tmp1));

    // Prepare the solution of the linear system (e.g. factorize)
    if (s.linsolF_.nfact(m->lin_nz, m->mem_linsolF)) casadi_error("'jac' factorization failed");

    return 0;
  } catch(int flag) { // recoverable error
    return flag;
  } catch(std::exception& e) { // non-recoverable error
    uerr() << "psetup failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesSimulator::calculate_jac(CvodesSimMemory* m, double t, N_Vector x, N_Vector xdot) const {
  // Input buffers
  std::fill_n(m->arg, enum_traits<DynIn>::n_enum + enum_traits<DynOut>::n_enum, nullptr);
  m->arg[DYN_T] = &t;  // t
  m->arg[DYN_X] = NV_DATA_S(x);  // x
  m->arg[DYN_U] = m->u;  // u
  m->arg[DYN_P] = m->p;  // p
  m->arg[enum_traits<DynIn>::n_enum + DYN_ODE] = NV_DATA_S(xdot);  // ode

  // Output buffers
  std::fill_n(m->res, enum_traits<DynIn>::n_enum * enum_traits<DynOut>::n_enum, nullptr);
  m->res[DYN_X + enum_traits<DynIn>::n_enum * DYN_ODE] = m->jac_x;  // jac:ode:x
  m->res[DYN_U + enum_traits<DynIn>::n_enum * DYN_ODE] = m->jac_u;  // jac:ode:u
  m->res[DYN_P + enum_traits<DynIn>::n_enum * DYN_ODE] = m->jac_p;  // jac:ode:p

#if 0
  // Update input cache
  bool changed = update_jac_in(m);

  // Quick return if inputs unchanged from last call
  if (!changed) return 0;
#endif

  // Evaluate
  if (calc_function(m, "jac")) return 1;

  // Successful return
  return 0;
}

int CvodesSimulator::lsetup(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
    booleantype *jcurPtr, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3) {
  try {
    auto m = to_mem(cv_mem->cv_lmem);
    //auto& s = m->self;

    // Current time
    double t = cv_mem->cv_tn;

    // Scaling factor before J
    double gamma = cv_mem->cv_gamma;

    // Call the preconditioner setup function (which sets up the linear solver)
    if (psetup(t, x, xdot, FALSE, jcurPtr,
               gamma, static_cast<void*>(m), vtemp1, vtemp2, vtemp3)) return 1;

    return 0;
  } catch(int flag) { // recoverable error
    return flag;
  } catch(std::exception& e) { // non-recoverable error
    uerr() << "lsetup failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesSimulator::lsolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
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
    if (psolve(t, x, xdot, b, b, gamma, delta,
               lr, static_cast<void*>(m), nullptr)) return 1;

    return 0;
  } catch(int flag) { // recoverable error
    return flag;
  } catch(std::exception& e) { // non-recoverable error
    uerr() << "lsolve failed: " << e.what() << std::endl;
    return -1;
  }
}

int CvodesSimulator::sens_rhs(int Ns, realtype t, N_Vector x, N_Vector xdot, N_Vector *fwd_x,
    N_Vector *fwd_xdot, void *user_data, N_Vector tmp1, N_Vector tmp2) {
  try {
    auto m = to_mem(user_data);
    auto& s = m->self;
    // Calculate Jacobian
    if (s.calculate_jac(m, t, x, xdot)) casadi_error("'jac' calculation failed");
    // Perform Jacobian-vector multiplications (can be parallelized if it matters)
    for (int k = 0; k < Ns; ++k) {
      // Consistency checks
      casadi_assert(NV_LENGTH_S(fwd_x[k]) == s.nx_, "Dimension mismatch");
      casadi_assert(NV_LENGTH_S(fwd_xdot[k]) == s.nx_, "Dimension mismatch");
      // Get pointers to data
      double *fwd_x_k = NV_DATA_S(fwd_x[k]), *fwd_xdot_k = NV_DATA_S(fwd_xdot[k]);
      // Clear fwd_xdot
      casadi_clear(fwd_xdot_k, s.nx_);
      // Matrix-vector multiplications
      casadi_mv(m->jac_x, s.jac_x_sp_, fwd_x_k, fwd_xdot_k, false);
      casadi_mv(m->jac_u, s.jac_u_sp_, m->fwd_u + s.nu_ * (s.ng_ - 1) * k, fwd_xdot_k, false);
      casadi_mv(m->jac_p, s.jac_p_sp_, m->fwd_p, fwd_xdot_k, false);
    }
    return 0;
  } catch(int flag) { // recoverable error
    return flag;
  } catch(std::exception& e) { // non-recoverable error
    uerr() << "CvodesSimulator::sens_rhs failed: " << e.what() << std::endl;
    return -1;
  }
}

CvodesSimMemory::CvodesSimMemory(const CvodesSimulator& s) : self(s) {
  this->mem = nullptr;

  // Reset checkpoints counter
  this->ncheck = 0;
}

CvodesSimMemory::~CvodesSimMemory() {
  if (this->mem_linsolF >= 0) self.linsolF_.release(this->mem_linsolF);
  if (this->mem) CVodeFree(&this->mem);
}

} // namespace casadi
