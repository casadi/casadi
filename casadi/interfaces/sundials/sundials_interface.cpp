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


#include "sundials_interface.hpp"

#include "casadi/core/casadi_misc.hpp"

INPUTSCHEME(IntegratorInput)
OUTPUTSCHEME(IntegratorOutput)

namespace casadi {

SundialsInterface::SundialsInterface(const std::string& name, const Function& dae,
    double t0, const std::vector<double>& tout)
    : Integrator(name, dae, t0, tout) {
}

SundialsInterface::~SundialsInterface() {
}

const Options SundialsInterface::options_
= {{&Integrator::options_},
    {{"max_num_steps",
      {OT_INT,
      "Maximum number of integrator steps"}},
    {"reltol",
      {OT_DOUBLE,
      "Relative tolerence for the IVP solution"}},
    {"abstol",
      {OT_DOUBLE,
      "Absolute tolerence for the IVP solution"}},
    {"newton_scheme",
      {OT_STRING,
      "Linear solver scheme in the Newton method: DIRECT|gmres|bcgstab|tfqmr"}},
    {"max_krylov",
      {OT_INT,
      "Maximum Krylov subspace size"}},
    {"sensitivity_method",
      {OT_STRING,
      "Sensitivity method: SIMULTANEOUS|staggered"}},
    {"max_multistep_order",
      {OT_INT,
      "Maximum order for the (variable-order) multistep method"}},
    {"use_preconditioner",
      {OT_BOOL,
      "Precondition the iterative solver [default: true]"}},
    {"stop_at_end",
      {OT_BOOL,
      "[DEPRECATED] Stop the integrator at the end of the interval"}},
    {"disable_internal_warnings",
      {OT_BOOL,
      "Disable SUNDIALS internal warning messages"}},
    {"quad_err_con",
      {OT_BOOL,
      "Should the quadratures affect the step size control"}},
    {"fsens_err_con",
      {OT_BOOL,
      "include the forward sensitivities in all error controls"}},
    {"steps_per_checkpoint",
      {OT_INT,
      "Number of steps between two consecutive checkpoints"}},
    {"interpolation_type",
      {OT_STRING,
      "Type of interpolation for the adjoint sensitivities"}},
    {"linear_solver",
      {OT_STRING,
      "A custom linear solver creator function [default: qr]"}},
    {"linear_solver_options",
      {OT_DICT,
      "Options to be passed to the linear solver"}},
    {"second_order_correction",
      {OT_BOOL,
      "Second order correction in the augmented system Jacobian [true]"}},
    {"step0",
      {OT_DOUBLE,
      "initial step size [default: 0/estimated]"}},
    {"max_step_size",
      {OT_DOUBLE,
      "Max step size [default: 0/inf]"}},
    {"max_order",
      {OT_DOUBLE,
      "Maximum order"}},
    {"nonlin_conv_coeff",
      {OT_DOUBLE,
      "Coefficient in the nonlinear convergence test"}},
    {"scale_abstol",
     {OT_BOOL,
      "Scale absolute tolerance by nominal value"}}
    }
};

void SundialsInterface::init(const Dict& opts) {
  // Call the base class method
  Integrator::init(opts);

  // Default options
  abstol_ = 1e-8;
  reltol_ = 1e-6;
  max_num_steps_ = 10000;
  stop_at_end_ = true;
  use_precon_ = true;
  max_krylov_ = 10;
  linear_solver_ = "qr";
  std::string newton_scheme = "direct";
  quad_err_con_ = false;
  std::string interpolation_type = "hermite";
  steps_per_checkpoint_ = 20;
  disable_internal_warnings_ = false;
  max_multistep_order_ = 5;
  second_order_correction_ = true;
  step0_ = 0;
  max_step_size_ = 0;
  max_order_ = 0;
  nonlin_conv_coeff_ = 0;
  scale_abstol_ = false;

  // Read options
  for (auto&& op : opts) {
    if (op.first=="abstol") {
      abstol_ = op.second;
    } else if (op.first=="reltol") {
      reltol_ = op.second;
    } else if (op.first=="max_num_steps") {
      max_num_steps_ = op.second;
    } else if (op.first=="stop_at_end") {
      stop_at_end_ = op.second;
      if (!stop_at_end_) {
        casadi_warning("The 'stop_at_end' option has been deprecated and is currently ignored");
      }
    } else if (op.first=="use_preconditioner") {
      use_precon_ = op.second;
    } else if (op.first=="max_krylov") {
      max_krylov_ = op.second;
    } else if (op.first=="newton_scheme") {
      newton_scheme = op.second.to_string();
    } else if (op.first=="linear_solver") {
      linear_solver_ = op.second.to_string();
    } else if (op.first=="linear_solver_options") {
      linear_solver_options_ = op.second;
    } else if (op.first=="quad_err_con") {
      quad_err_con_ = op.second;
    } else if (op.first=="interpolation_type") {
      interpolation_type = op.second.to_string();
    } else if (op.first=="steps_per_checkpoint") {
      steps_per_checkpoint_ = op.second;
    } else if (op.first=="disable_internal_warnings") {
      disable_internal_warnings_ = op.second;
    } else if (op.first=="max_multistep_order") {
      max_multistep_order_ = op.second;
    } else if (op.first=="second_order_correction") {
      second_order_correction_ = op.second;
    } else if (op.first=="step0") {
      step0_ = op.second;
    } else if (op.first=="max_step_size") {
      max_step_size_ = op.second;
    } else if (op.first=="max_order") {
      max_order_ = op.second;
    } else if (op.first=="nonlin_conv_coeff") {
      nonlin_conv_coeff_ = op.second;
    } else if (op.first=="scale_abstol") {
      scale_abstol_ = op.second;
    }
  }

  // Type of Newton scheme
  if (newton_scheme=="direct") {
    newton_scheme_ = SD_DIRECT;
  } else if (newton_scheme=="gmres") {
    newton_scheme_ = SD_GMRES;
  } else if (newton_scheme=="bcgstab") {
    newton_scheme_ = SD_BCGSTAB;
  } else if (newton_scheme=="tfqmr") {
    newton_scheme_ = SD_TFQMR;
  } else {
    casadi_error("Unknown Newton scheme: " + newton_scheme);
  }

  // Interpolation_type
  if (interpolation_type=="hermite") {
    interp_ = SD_HERMITE;
  } else if (interpolation_type=="polynomial") {
    interp_ = SD_POLYNOMIAL;
  } else {
    casadi_error("Unknown interpolation type: " + interpolation_type);
  }

  // If derivative, use Jacobian from non-augmented system if possible
  SundialsInterface* d = 0;
  if (nfwd_ > 0 && !derivative_of_.is_null()) {
    d = derivative_of_.get<SundialsInterface>();
    casadi_assert_dev(d != nullptr);
  }

  // Get Jacobian function, forward problem
  Function jacF;
  Sparsity jacF_sp;
  if (d == 0) {
    // New Jacobian function
    jacF = create_function("jacF", {"t", "x", "z", "p", "u"},
      {"jac:ode:x", "jac:alg:x", "jac:ode:z", "jac:alg:z"});
    jacF_sp = jacF.sparsity_out(JACF_ODE_X) + Sparsity::diag(nx1_);
    if (nz_ > 0) {
      jacF_sp = horzcat(vertcat(jacF_sp, jacF.sparsity_out(JACF_ALG_X)),
        vertcat(jacF.sparsity_out(JACF_ODE_Z), jacF.sparsity_out(JACF_ALG_Z)));
    }
  } else {
    // Reuse existing Jacobian function
    jacF = d->get_function("jacF");
    set_function(jacF, jacF.name(), true);
    linsolF_ = d->linsolF_;
    jacF_sp = linsolF_.sparsity();
  }
  alloc_w(jacF_sp.nnz(), true);  // jacF

  // Linear solver for forward problem
  if (linsolF_.is_null()) {
    linsolF_ = Linsol("linsolF", linear_solver_, jacF_sp, linear_solver_options_);
  }

  // Attach functions to calculate DAE and quadrature RHS all-at-once
  if (nfwd_ > 0) {
    create_forward("daeF", nfwd_);
    if (nq_ > 0) create_forward("quadF", nfwd_);
    if (nadj_ > 0) {
      create_forward("daeB", nfwd_);
      if (nrq_ > 0 || nuq_ > 0) create_forward("quadB", nfwd_);
    }
  }

  // Attach functions for jacobian information, foward problem
  if (newton_scheme_ != SD_DIRECT || (nfwd_ > 0 && second_order_correction_)) {
    create_function("jtimesF", {"t", "x", "z", "p", "u", "fwd:x", "fwd:z"},
      {"fwd:ode", "fwd:alg"});
    if (nfwd_ > 0) {
      create_forward("jtimesF", nfwd_);
    }
  }


  // For Jacobian calculation
  alloc_w(jacF.nnz_out(JACF_ODE_X), true);  // jac_ode_x
  alloc_w(jacF.nnz_out(JACF_ALG_X), true);  // jac_alg_x
  alloc_w(jacF.nnz_out(JACF_ODE_Z), true);  // jac_ode_z
  alloc_w(jacF.nnz_out(JACF_ALG_Z), true);  // jac_alg_z

  // Transposing the Jacobian (for calculating jacB)
  // This will be unnecessary once linsolF_ is used for both forward and adjoint
  // cf. #3047
  alloc_w(nx_ + nz_); // casadi_trans
  alloc_iw(nx_ + nz_); // casadi_trans
}

void SundialsInterface::set_work(void* mem, const double**& arg, double**& res,
    casadi_int*& iw, double*& w) const {
  auto m = static_cast<SundialsMemory*>(mem);

  // Set work in base classes
  Integrator::set_work(mem, arg, res, iw, w);

  // Work vectors
  m->jacF = w; w += linsolF_.sparsity().nnz();

  // Work vectors
  const Function& jacF = get_function("jacF");
  m->jac_ode_x = w; w += jacF.nnz_out(JACF_ODE_X);
  m->jac_alg_x = w; w += jacF.nnz_out(JACF_ALG_X);
  m->jac_ode_z = w; w += jacF.nnz_out(JACF_ODE_Z);
  m->jac_alg_z = w; w += jacF.nnz_out(JACF_ALG_Z);
}

int SundialsInterface::init_mem(void* mem) const {
  if (Integrator::init_mem(mem)) return 1;
  auto m = static_cast<SundialsMemory*>(mem);

  // Allocate NVectors
  m->xz = N_VNew_Serial(nx_ + nz_);
  m->q = N_VNew_Serial(nq_);
  m->rxz = N_VNew_Serial(nrx_ + nrz_);
  m->ruq = N_VNew_Serial(nrq_ + nuq_);

  // Absolute tolerances as NVector
  if (scale_abstol_) {
    // Allocate NVector
    m->abstolv = N_VNew_Serial(nx_ + nz_);
    // Get pointer to data
    double* abstolv = NV_DATA_S(m->abstolv);
    // States
    for (casadi_int d = 0; d <= nfwd_; ++d) {
      for (casadi_int i = 0; i < nx1_; ++i) *abstolv++ = abstol_ * nom_x_[i];
    }
    // Algebraic variables
    for (casadi_int d = 0; d <= nfwd_; ++d) {
      for (casadi_int i = 0; i < nz1_; ++i) *abstolv++ = abstol_ * nom_z_[i];
    }
    // Consistency check
    casadi_assert_dev(abstolv == NV_DATA_S(m->abstolv) + nx_ + nz_);
  } else {
    m->abstolv = nullptr;
  }

  m->mem_linsolF = linsolF_.checkout();

  // Reset stats
  reset_stats(m);

  return 0;
}

void SundialsInterface::reset(IntegratorMemory* mem) const {
  auto m = static_cast<SundialsMemory*>(mem);

  // Reset the base classes
  Integrator::reset(mem);

  // Reset stats
  reset_stats(m);

  // Set parameters
  casadi_copy(p, np_, m->p);

  // Set controls
  casadi_copy(u, nu_, m->u);

  // Set the state
  casadi_copy(m->x, nx_, NV_DATA_S(m->xz));
  casadi_copy(m->z, nz_, NV_DATA_S(m->xz) + nx_);

  // Reset summation states
  N_VConst(0., m->q);
}

void SundialsInterface::reset_stats(SundialsMemory* m) const {
  // Reset stats, forward problem
  m->nsteps = m->nfevals = m->nlinsetups = m->netfails = 0;
  m->qlast = m->qcur = -1;
  m->hinused = m->hlast = m->hcur = m->tcur = casadi::nan;
  m->nniters = m->nncfails = 0;

  // Reset stats, backward problem
  m->nstepsB = m->nfevalsB = m->nlinsetupsB = m->netfailsB = 0;
  m->qlastB = m->qcurB = -1;
  m->hinusedB = m->hlastB = m->hcurB = m->tcurB = casadi::nan;
  m->nnitersB = m->nncfailsB = 0;

  // Set offsets to zero
  save_offsets(m);
}

void SundialsInterface::save_offsets(SundialsMemory* m) const {
  // Retrieve stats offset, backward problem
  m->nstepsB_off = m->nstepsB;
  m->nfevalsB_off = m->nfevalsB;
  m->nlinsetupsB_off = m->nlinsetupsB;
  m->netfailsB_off = m->netfailsB;
  m->nnitersB_off = m->nnitersB;
  m->nncfailsB_off = m->nncfailsB;
}

void SundialsInterface::add_offsets(SundialsMemory* m) const {
  // Add stats offsets, backward problem
  m->nstepsB += m->nstepsB_off;
  m->nfevalsB += m->nfevalsB_off;
  m->nlinsetupsB += m->nlinsetupsB_off;
  m->netfailsB += m->netfailsB_off;
  m->nnitersB += m->nnitersB_off;
  m->nncfailsB += m->nncfailsB_off;

}

void SundialsInterface::resetB(IntegratorMemory* mem) const {
  auto m = static_cast<SundialsMemory*>(mem);

  // Clear seeds
  casadi_clear(m->adj_q, nrp_);
  casadi_clear(NV_DATA_S(m->rxz), nrx_ + nrz_);

  // Reset summation states
  N_VConst(0., m->ruq);
}

void SundialsInterface::impulseB(IntegratorMemory* mem,
    const double* adj_x, const double* adj_z, const double* adj_q) const {
  auto m = static_cast<SundialsMemory*>(mem);

  // Add impulse to backward parameters
  casadi_axpy(nrp_, 1., adj_q, m->adj_q);

  // Add impulse to backward state
  casadi_axpy(nrx_, 1., adj_x, NV_DATA_S(m->rxz));

  // Add impulse to algebraic variables:
  // If nonzero, this has to be propagated to an impulse in backward state
  // casadi_copy(adj_z, nrz_, NV_DATA_S(m->rxz) + nrx_);
  casadi_axpy(nrz_, 1., adj_z, NV_DATA_S(m->rxz) + nrx_);
}

SundialsMemory::SundialsMemory() {
  this->xz  = nullptr;
  this->q = nullptr;
  this->rxz = nullptr;
  this->ruq = nullptr;
  this->first_callB = true;
  this->abstolv = nullptr;
  this->mem_linsolF = -1;
}

SundialsMemory::~SundialsMemory() {
  if (this->xz) N_VDestroy_Serial(this->xz);
  if (this->q) N_VDestroy_Serial(this->q);
  if (this->rxz) N_VDestroy_Serial(this->rxz);
  if (this->ruq) N_VDestroy_Serial(this->ruq);
  if (this->abstolv) N_VDestroy_Serial(this->abstolv);
}

Dict SundialsInterface::get_stats(void* mem) const {
  Dict stats = Integrator::get_stats(mem);
  auto m = static_cast<SundialsMemory*>(mem);

  // Counters, forward problem
  stats["nsteps"] = static_cast<casadi_int>(m->nsteps);
  stats["nfevals"] = static_cast<casadi_int>(m->nfevals);
  stats["nlinsetups"] = static_cast<casadi_int>(m->nlinsetups);
  stats["netfails"] = static_cast<casadi_int>(m->netfails);
  stats["qlast"] = m->qlast;
  stats["qcur"] = m->qcur;
  stats["hinused"] = m->hinused;
  stats["hlast"] = m->hlast;
  stats["hcur"] = m->hcur;
  stats["tcur"] = m->tcur;
  stats["nniters"] = static_cast<casadi_int>(m->nniters);
  stats["nncfails"] = static_cast<casadi_int>(m->nncfails);

  // Counters, backward problem
  stats["nstepsB"] = static_cast<casadi_int>(m->nstepsB);
  stats["nfevalsB"] = static_cast<casadi_int>(m->nfevalsB);
  stats["nlinsetupsB"] = static_cast<casadi_int>(m->nlinsetupsB);
  stats["netfailsB"] = static_cast<casadi_int>(m->netfailsB);
  stats["qlastB"] = m->qlastB;
  stats["qcurB"] = m->qcurB;
  stats["hinusedB"] = m->hinusedB;
  stats["hlastB"] = m->hlastB;
  stats["hcurB"] = m->hcurB;
  stats["tcurB"] = m->tcurB;
  stats["nnitersB"] = static_cast<casadi_int>(m->nnitersB);
  stats["nncfailsB"] = static_cast<casadi_int>(m->nncfailsB);
  return stats;
}

void SundialsInterface::print_stats(IntegratorMemory* mem) const {
  auto m = to_mem(mem);
  print("FORWARD INTEGRATION:\n");
  print("Number of steps taken by SUNDIALS: %ld\n", m->nsteps);
  print("Number of calls to the user's f function: %ld\n", m->nfevals);
  print("Number of calls made to the linear solver setup function: %ld\n", m->nlinsetups);
  print("Number of error test failures: %ld\n", m->netfails);
  print("Method order used on the last internal step: %d\n", m->qlast);
  print("Method order to be used on the next internal step: %d\n", m->qcur);
  print("Actual value of initial step size: %g\n", m->hinused);
  print("Step size taken on the last internal step: %g\n", m->hlast);
  print("Step size to be attempted on the next internal step: %g\n", m->hcur);
  print("Current internal time reached: %g\n", m->tcur);
  print("Number of nonlinear iterations performed: %ld\n", m->nniters);
  print("Number of nonlinear convergence failures: %ld\n", m->nncfails);
  if (nrx_>0) {
    print("BACKWARD INTEGRATION:\n");
    print("Number of steps taken by SUNDIALS: %ld\n", m->nstepsB);
    print("Number of calls to the user's f function: %ld\n", m->nfevalsB);
    print("Number of calls made to the linear solver setup function: %ld\n", m->nlinsetupsB);
    print("Number of error test failures: %ld\n", m->netfailsB);
    print("Method order used on the last internal step: %d\n" , m->qlastB);
    print("Method order to be used on the next internal step: %d\n", m->qcurB);
    print("Actual value of initial step size: %g\n", m->hinusedB);
    print("Step size taken on the last internal step: %g\n", m->hlastB);
    print("Step size to be attempted on the next internal step: %g\n", m->hcurB);
    print("Current internal time reached: %g\n", m->tcurB);
    print("Number of nonlinear iterations performed: %ld\n", m->nnitersB);
    print("Number of nonlinear convergence failures: %ld\n", m->nncfailsB);
  }
  print("\n");
}

SundialsInterface::SundialsInterface(DeserializingStream& s) : Integrator(s) {
  int version = s.version("SundialsInterface", 1, 2);
  s.unpack("SundialsInterface::abstol", abstol_);
  s.unpack("SundialsInterface::reltol", reltol_);
  s.unpack("SundialsInterface::max_num_steps", max_num_steps_);
  s.unpack("SundialsInterface::stop_at_end", stop_at_end_);
  s.unpack("SundialsInterface::quad_err_con", quad_err_con_);
  s.unpack("SundialsInterface::steps_per_checkpoint", steps_per_checkpoint_);
  s.unpack("SundialsInterface::disable_internal_warnings", disable_internal_warnings_);
  s.unpack("SundialsInterface::max_multistep_order", max_multistep_order_);
  s.unpack("SundialsInterface::linear_solver", linear_solver_);
  s.unpack("SundialsInterface::linear_solver_options", linear_solver_options_);

  s.unpack("SundialsInterface::max_krylov", max_krylov_);
  s.unpack("SundialsInterface::use_precon", use_precon_);
  s.unpack("SundialsInterface::second_order_correction", second_order_correction_);

  s.unpack("SundialsInterface::step0", step0_);
  if (version>=2) {
    s.unpack("SundialsInterface::max_step_size", max_step_size_);
  } else {
    max_step_size_ = 0;
  }

  s.unpack("SundialsInterface::nonlin_conv_coeff", nonlin_conv_coeff_);
  s.unpack("SundialsInterface::max_order", max_order_);
  s.unpack("SundialsInterface::scale_abstol", scale_abstol_);

  s.unpack("SundialsInterface::linsolF", linsolF_);

  int newton_scheme;
  s.unpack("SundialsInterface::newton_scheme", newton_scheme);
  newton_scheme_ = static_cast<NewtonScheme>(newton_scheme);

  int interp;
  s.unpack("SundialsInterface::interp", interp);
  interp_ = static_cast<InterpType>(interp);

}

void SundialsInterface::serialize_body(SerializingStream &s) const {
  Integrator::serialize_body(s);
  s.version("SundialsInterface", 2);
  s.pack("SundialsInterface::abstol", abstol_);
  s.pack("SundialsInterface::reltol", reltol_);
  s.pack("SundialsInterface::max_num_steps", max_num_steps_);
  s.pack("SundialsInterface::stop_at_end", stop_at_end_);
  s.pack("SundialsInterface::quad_err_con", quad_err_con_);
  s.pack("SundialsInterface::steps_per_checkpoint", steps_per_checkpoint_);
  s.pack("SundialsInterface::disable_internal_warnings", disable_internal_warnings_);
  s.pack("SundialsInterface::max_multistep_order", max_multistep_order_);

  s.pack("SundialsInterface::linear_solver", linear_solver_);
  s.pack("SundialsInterface::linear_solver_options", linear_solver_options_);
  s.pack("SundialsInterface::max_krylov", max_krylov_);
  s.pack("SundialsInterface::use_precon", use_precon_);
  s.pack("SundialsInterface::second_order_correction", second_order_correction_);

  s.pack("SundialsInterface::step0", step0_);
  s.pack("SundialsInterface::max_step_size", max_step_size_);

  s.pack("SundialsInterface::nonlin_conv_coeff", nonlin_conv_coeff_);
  s.pack("SundialsInterface::max_order", max_order_);
  s.pack("SundialsInterface::scale_abstol", scale_abstol_);

  s.pack("SundialsInterface::linsolF", linsolF_);

  s.pack("SundialsInterface::newton_scheme", static_cast<int>(newton_scheme_));
  s.pack("SundialsInterface::interp", static_cast<int>(interp_));
}

int SundialsInterface::calc_daeF(SundialsMemory* m, double t, const double* x, const double* z,
    double* ode, double* alg) const {
  // Evaluate nondifferentiated
  m->arg[DYN_T] = &t;  // t
  m->arg[DYN_X] = x;  // x
  m->arg[DYN_Z] = z;  // z
  m->arg[DYN_P] = m->p;  // p
  m->arg[DYN_U] = m->u;  // u
  m->res[DAE_ODE] = ode;  // ode
  m->res[DAE_ALG] = alg;  // alg
  if (calc_function(m, "daeF")) return 1;
  // Evaluate sensitivities
  if (nfwd_ > 0) {
    m->arg[DYN_NUM_IN + DAE_ODE] = ode;  // out:ode
    m->arg[DYN_NUM_IN + DAE_ALG] = alg;  // out:alg
    m->arg[DYN_NUM_IN + DAE_NUM_OUT + DYN_T] = 0;  // fwd:t
    m->arg[DYN_NUM_IN + DAE_NUM_OUT + DYN_X] = x + nx1_;  // fwd:x
    m->arg[DYN_NUM_IN + DAE_NUM_OUT + DYN_Z] = z ? z + nz1_ : 0;  // fwd:z
    m->arg[DYN_NUM_IN + DAE_NUM_OUT + DYN_P] = m->p + np1_;  // fwd:p
    m->arg[DYN_NUM_IN + DAE_NUM_OUT + DYN_U] = m->u + nu1_;  // fwd:u
    m->res[DAE_ODE] = ode ? ode + nx1_ : 0;  // fwd:ode
    m->res[DAE_ALG] = alg ? alg + nz1_ : 0;  // fwd:alg
    if (calc_function(m, forward_name("daeF", nfwd_))) return 1;
  }
  return 0;
}

int SundialsInterface::calc_daeB(SundialsMemory* m, double t, const double* x, const double* z,
    const double* adj_ode, const double* adj_alg, const double* adj_quad,
    double* adj_x, double* adj_z) const {
  // Evaluate nondifferentiated
  m->arg[BDYN_T] = &t;  // t
  m->arg[BDYN_X] = x;  // x
  m->arg[BDYN_Z] = z;  // z
  m->arg[BDYN_P] = m->p;  // p
  m->arg[BDYN_U] = m->u;  // u
  m->arg[BDYN_OUT_ODE] = nullptr;  // out_ode
  m->arg[BDYN_OUT_ALG] = nullptr;  // out_alg
  m->arg[BDYN_OUT_QUAD] = nullptr;  // out_quad
  m->arg[BDYN_OUT_ZERO] = nullptr;  // out_zero
  m->arg[BDYN_ADJ_ODE] = adj_ode;  // adj_ode
  m->arg[BDYN_ADJ_ALG] = adj_alg;  // adj_alg
  m->arg[BDYN_ADJ_QUAD] = adj_quad;  // adj_quad
  m->arg[BDYN_ADJ_ZERO] = nullptr;  // adj_zero
  m->res[BDAE_ADJ_X] = adj_x;  // adj_x
  m->res[BDAE_ADJ_Z] = adj_z;  // adj_z
  if (calc_function(m, "daeB")) return 1;
  // Evaluate sensitivities
  if (nfwd_ > 0) {
    m->arg[BDYN_NUM_IN + BDAE_ADJ_X] = adj_x;  // out:adj_x
    m->arg[BDYN_NUM_IN + BDAE_ADJ_Z] = adj_z;  // out:adj_z
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_T] = 0;  // fwd:t
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_X] = x ? x + nx1_ : x;  // fwd:x
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_Z] = z ? z + nz1_ : z;  // fwd:z
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_P] = m->p + np1_;  // fwd:p
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_U] = m->u + nu1_;  // fwd:u
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_OUT_ODE] = nullptr;  // fwd:out_ode
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_OUT_ALG] = nullptr;  // fwd:out_alg
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_OUT_QUAD] = nullptr;  // fwd:out_quad
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_OUT_ZERO] = nullptr;  // fwd:out_zero
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_ADJ_ODE] =
      adj_ode ? adj_ode + nrx1_ * nadj_ : 0;  // fwd:adj_ode
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_ADJ_ALG] =
      adj_alg ? adj_alg + nrz1_ * nadj_ : 0;  // fwd:adj_alg
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_ADJ_QUAD] =
      adj_quad ? adj_quad + nrp1_ * nadj_ : 0;  // fwd:adj_quad
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_ADJ_ZERO] = nullptr;  // fwd:adj_zero
    m->res[BDAE_ADJ_X] = adj_x ? adj_x + nrx1_ * nadj_ : 0;  // fwd:adj_x
    m->res[BDAE_ADJ_Z] = adj_z ? adj_z + nrz1_ * nadj_ : 0;  // fwd:adj_z
    if (calc_function(m, forward_name("daeB", nfwd_))) return 1;
  }
  return 0;
}

int SundialsInterface::calc_quadF(SundialsMemory* m, double t, const double* x, const double* z,
    double* quad) const {
  m->arg[DYN_T] = &t;  // t
  m->arg[DYN_X] = x;  // x
  m->arg[DYN_Z] = z;  // z
  m->arg[DYN_P] = m->p;  // p
  m->arg[DYN_U] = m->u;  // u
  m->res[QUAD_QUAD] = quad;  // quad
  if (calc_function(m, "quadF")) return 1;
  // Evaluate sensitivities
  if (nfwd_ > 0) {
    m->arg[DYN_NUM_IN + QUAD_QUAD] = quad;  // out:quad
    m->arg[DYN_NUM_IN + QUAD_NUM_OUT + DYN_T] = 0;  // fwd:t
    m->arg[DYN_NUM_IN + QUAD_NUM_OUT + DYN_X] = x + nx1_;  // fwd:x
    m->arg[DYN_NUM_IN + QUAD_NUM_OUT + DYN_Z] = z ? z + nz1_ : 0;  // fwd:z
    m->arg[DYN_NUM_IN + QUAD_NUM_OUT + DYN_P] = m->p + np1_;  // fwd:p
    m->arg[DYN_NUM_IN + QUAD_NUM_OUT + DYN_U] = m->u + nu1_;  // fwd:u
    m->res[QUAD_QUAD] = quad ? quad + nq1_ : 0;  // fwd:quad
    if (calc_function(m, forward_name("quadF", nfwd_))) return 1;
  }
  return 0;
}

int SundialsInterface::calc_quadB(SundialsMemory* m, double t, const double* x, const double* z,
    const double* adj_ode, const double* adj_alg, double* adj_p, double* adj_u) const {
  // Evaluate nondifferentiated
  m->arg[BDYN_T] = &t;  // t
  m->arg[BDYN_X] = x;  // x
  m->arg[BDYN_Z] = z;  // z
  m->arg[BDYN_P] = m->p;  // p
  m->arg[BDYN_U] = m->u;  // u
  m->arg[BDYN_OUT_ODE] = nullptr;  // out_ode
  m->arg[BDYN_OUT_ALG] = nullptr;  // out_alg
  m->arg[BDYN_OUT_QUAD] = nullptr;  // out_quad
  m->arg[BDYN_OUT_ZERO] = nullptr;  // out_zero
  m->arg[BDYN_ADJ_ODE] = adj_ode;  // adj_ode
  m->arg[BDYN_ADJ_ALG] = adj_alg;  // adj_alg
  m->arg[BDYN_ADJ_QUAD] = m->adj_q;  // adj_quad
  m->arg[BDYN_ADJ_ZERO] = nullptr;  // adj_zero
  m->res[BQUAD_ADJ_P] = adj_p;  // adj_p
  m->res[BQUAD_ADJ_U] = adj_u;  // adj_u
  if (calc_function(m, "quadB")) return 1;
  // Evaluate sensitivities
  if (nfwd_ > 0) {
    m->arg[BDYN_NUM_IN + BQUAD_ADJ_P] = adj_p;  // out:adj_p
    m->arg[BDYN_NUM_IN + BQUAD_ADJ_U] = adj_u;  // out:adj_u
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_T] = 0;  // fwd:t
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_X] = x ? x + nx1_ : 0;  // fwd:x
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_Z] = z ? z + nz1_ : 0;  // fwd:z
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_P] = m->p + np1_;  // fwd:p
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_U] = m->u + nu1_;  // fwd:u
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_OUT_ODE] = nullptr;  // fwd:out_ode
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_OUT_ALG] = nullptr;  // fwd:out_alg
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_OUT_QUAD] = nullptr;  // fwd:out_quad
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_OUT_ZERO] = nullptr;  // fwd:out_zero
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_ADJ_ODE] =
      adj_ode ? adj_ode + nrx1_ * nadj_ : 0;  // fwd:adj_ode
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_ADJ_ALG] =
      adj_alg ? adj_alg + nrz1_ * nadj_ : 0;  // fwd:adj_alg
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_ADJ_QUAD] = m->adj_q + nrp1_ * nadj_;  // fwd:adj_quad
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_ADJ_ZERO] = nullptr;  // fwd:adj_zero
    m->res[BQUAD_ADJ_P] = adj_p + nrq1_ * nadj_;  // fwd:adj_p
    m->res[BQUAD_ADJ_U] = adj_u + nuq1_ * nadj_;  // fwd:adj_u
    if (calc_function(m, forward_name("quadB", nfwd_))) return 1;
  }
  return 0;
}

int SundialsInterface::calc_jtimesF(SundialsMemory* m, double t, const double* x, const double* z,
    const double* fwd_x, const double* fwd_z, double* fwd_ode, double* fwd_alg) const {
  // Evaluate nondifferentiated
  m->arg[JTIMESF_T] = &t;  // t
  m->arg[JTIMESF_X] = x;  // x
  m->arg[JTIMESF_Z] = z;  // z
  m->arg[JTIMESF_P] = m->p;  // p
  m->arg[JTIMESF_U] = m->u;  // u
  m->arg[JTIMESF_FWD_X] = fwd_x;  // fwd:x
  m->arg[JTIMESF_FWD_Z] = fwd_z;  // fwd:z
  m->res[JTIMESF_FWD_ODE] = fwd_ode;  // fwd:ode
  m->res[JTIMESF_FWD_ALG] = fwd_alg;  // fwd:alg
  if (calc_function(m, "jtimesF")) return 1;
  // Evaluate sensitivities
  if (nfwd_ > 0) {
    m->arg[JTIMESF_NUM_IN + JTIMESF_FWD_ODE] = fwd_ode;  // out:fwd:ode
    m->arg[JTIMESF_NUM_IN + JTIMESF_FWD_ALG] = fwd_alg;  // out:fwd:alg
    m->arg[JTIMESF_NUM_IN + JTIMESF_NUM_OUT + JTIMESF_T] = 0;  // fwd:t
    m->arg[JTIMESF_NUM_IN + JTIMESF_NUM_OUT + JTIMESF_X] = x + nx1_;  // fwd:x
    m->arg[JTIMESF_NUM_IN + JTIMESF_NUM_OUT + JTIMESF_Z] = z + nz1_;  // fwd:z
    m->arg[JTIMESF_NUM_IN + JTIMESF_NUM_OUT + JTIMESF_P] = m->p + np1_;  // fwd:p
    m->arg[JTIMESF_NUM_IN + JTIMESF_NUM_OUT + JTIMESF_U] = m->u + nu1_;  // fwd:u
    m->arg[JTIMESF_NUM_IN + JTIMESF_NUM_OUT + JTIMESF_FWD_X] = fwd_x + nx1_;  // fwd:fwd:x
    m->arg[JTIMESF_NUM_IN + JTIMESF_NUM_OUT + JTIMESF_FWD_Z] = fwd_z + nz1_;  // fwd:fwd:z
    m->res[JTIMESF_FWD_ODE] = fwd_ode + nx1_;  // fwd:fwd:ode
    m->res[JTIMESF_FWD_ALG] = fwd_alg + nz1_;  // fwd:fwd:alg
    if (calc_function(m, forward_name("jtimesF", nfwd_))) return 1;
  }
  // Successful return
  return 0;
}

int SundialsInterface::calc_jacF(SundialsMemory* m, double t, const double* x, const double* z,
    double* jac_ode_x, double* jac_alg_x, double* jac_ode_z, double* jac_alg_z) const {
  // Calculate Jacobian
  m->arg[DYN_T] = &t;
  m->arg[DYN_X] = x;
  m->arg[DYN_Z] = z;
  m->arg[DYN_P] = m->p;
  m->arg[DYN_U] = m->u;
  m->res[JACF_ODE_X] = jac_ode_x;
  m->res[JACF_ALG_X] = jac_alg_x;
  m->res[JACF_ODE_Z] = jac_ode_z;
  m->res[JACF_ALG_Z] = jac_alg_z;
  return calc_function(m, "jacF");
}

} // namespace casadi
