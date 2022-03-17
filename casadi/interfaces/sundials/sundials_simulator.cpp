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


#include "sundials_simulator.hpp"

#include "casadi/core/casadi_misc.hpp"

INPUTSCHEME(SimulatorInput)
OUTPUTSCHEME(SimulatorOutput)

namespace casadi {

std::string to_string(NewtonScheme v) {
  switch (v) {
  case NewtonScheme::DIRECT: return "direct";
  case NewtonScheme::GMRES: return "gmres";
  case NewtonScheme::BCGSTAB: return "bcgstab";
  case NewtonScheme::TFQMR: return "tfqmr";
  default: break;
  }
  return "";
}

std::string to_string(InterpType v) {
  switch (v) {
  case InterpType::POLYNOMIAL: return "polynomial";
  case InterpType::HERMITE: return "hermite";
  default: break;
  }
  return "";
}

SundialsSimulator::SundialsSimulator(const std::string& name, const Function& dae,
    const std::vector<double>& grid)
  : Simulator(name, dae, grid) {
}

SundialsSimulator::~SundialsSimulator() {
}

const Options SundialsSimulator::options_
= {{&Simulator::options_},
   {{"max_num_steps",
     {OT_INT,
      "Maximum number of simulator steps"}},
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
      "Stop the simulator at the end of the interval"}},
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

void SundialsSimulator::init(const Dict& opts) {
  // Call the base class method
  Simulator::init(opts);

  // Default options
  abstol_ = 1e-5;
  reltol_ = 1e-5;
  max_num_steps_ = 5000;
  stop_at_end_ = true;
  use_precon_ = true;
  max_krylov_ = 10;
  linear_solver_ = "qr";
  newton_scheme_ = NewtonScheme::DIRECT;
  quad_err_con_ = false;
  interp_ = InterpType::HERMITE;
  steps_per_checkpoint_ = 20;
  disable_internal_warnings_ = false;
  max_multistep_order_ = 5;
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
    } else if (op.first=="use_preconditioner") {
      use_precon_ = op.second;
    } else if (op.first=="max_krylov") {
      max_krylov_ = op.second;
    } else if (op.first=="newton_scheme") {
      newton_scheme_ = to_enum<NewtonScheme>(op.second.to_string());
    } else if (op.first=="linear_solver") {
      linear_solver_ = op.second.to_string();
    } else if (op.first=="linear_solver_options") {
      linear_solver_options_ = op.second;
    } else if (op.first=="quad_err_con") {
      quad_err_con_ = op.second;
    } else if (op.first=="interpolation_type") {
      interp_ = to_enum<InterpType>(op.second.to_string());
    } else if (op.first=="steps_per_checkpoint") {
      steps_per_checkpoint_ = op.second;
    } else if (op.first=="disable_internal_warnings") {
      disable_internal_warnings_ = op.second;
    } else if (op.first=="max_multistep_order") {
      max_multistep_order_ = op.second;
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

  // Get or create Jacobians and linear system solvers
  Function J = oracle_.jacobian();
  set_function(J, "jac", true);

  // Allocate work vectors
  // alloc_w(nx_ + nz_, true);  // xz
  alloc_w(2 * (nx_+nz_), true); // v1, v2

  // Allocate linear solvers
  Sparsity Jsp = J.sparsity_out("jac_ode_x");
  Jsp = Jsp + Sparsity::diag(Jsp.size());
  linsolF_ = Linsol("linsolF", linear_solver_, Jsp, linear_solver_options_);
  alloc_w(linsolF_.sparsity().nnz(), true);
}

int SundialsSimulator::init_mem(void* mem) const {
  if (Simulator::init_mem(mem)) return 1;
  auto m = static_cast<SundialsSimMemory*>(mem);

  // Allocate n-vectors
  m->xz = N_VNew_Serial(nx_+nz_);
  if (scale_abstol_) {
    m->abstolv = N_VNew_Serial(nx_+nz_);
  }

  m->mem_linsolF = linsolF_.checkout();

  return 0;
}

void SundialsSimulator::reset(SimulatorMemory* mem) const {
  auto m = static_cast<SundialsSimMemory*>(mem);
  // Set the state
  casadi_copy(m->xk, nx_, NV_DATA_S(m->xz));
  casadi_copy(m->zk, nz_, NV_DATA_S(m->xz) + nx_);
}

SundialsSimMemory::SundialsSimMemory() {
  // Set pointers to null
  this->xz  = nullptr;
  this->jac = nullptr;
  this->v1 = this->v2 = nullptr;
  this->abstolv  = nullptr;
  // Reset stats
  this->nsteps = 0;
  this->nfevals = 0;
  this->nlinsetups = 0;
  this->netfails = 0;
  this->qlast = 0;
  this->qcur = 0;
  this->hinused = 0;
  this->hlast = 0;
  this->hcur = 0;
  this->tcur = 0;
  this->nniters = 0;
  this->nncfails = 0;
  // Number of checkpoints stored so far
  this->ncheck = 0;
  // Unused
  this->mem_linsolF = -1;
}

SundialsSimMemory::~SundialsSimMemory() {
  if (this->xz) N_VDestroy_Serial(this->xz);
  if (this->abstolv) N_VDestroy_Serial(this->abstolv);
}

Dict SundialsSimulator::get_stats(void* mem) const {
  Dict stats = Simulator::get_stats(mem);
  auto m = static_cast<SundialsSimMemory*>(mem);

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
  return stats;
}

void SundialsSimulator::print_stats(SimulatorMemory* mem) const {
  auto m = to_mem(mem);
  print("Number of steps taken by SUNDIALS: %ld\n", m->nsteps);
  print("Number of calls to the user’s f function: %ld\n", m->nfevals);
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
  print("\n");
}

void SundialsSimulator::set_work(void* mem, const double**& arg, double**& res,
                              casadi_int*& iw, double*& w) const {
  auto m = static_cast<SundialsSimMemory*>(mem);

  // Set work in base classes
  Simulator::set_work(mem, arg, res, iw, w);

  // Work vectors
  m->v1 = w; w += nx_ + nz_;
  m->v2 = w; w += nx_ + nz_;
  m->jac = w; w += linsolF_.sparsity().nnz();
}

void SundialsSimulator::add_diag(const casadi_int* sp, double *nz, double v,
    const casadi_int* sp_new, double *w) {
  // Local variables
  casadi_int nrow, ncol, r, c, k;
  const casadi_int *colind, *row, *colind_new, *row_new;
  // Extract sparsities
  nrow = sp[0];
  ncol = sp[1];
  colind = sp + 2;
  row = colind + ncol + 1;
  colind_new = sp_new + 2;
  row_new = colind_new + ncol + 1;
  // Clear work vector
  for (r = 0; r < nrow; ++r) w[r] = 0;
  // Loop over columns in reverse order
  for (c = ncol; c-- > 0; ) {
    // Copy content of column
    for (k = colind[c]; k < colind[c + 1]; ++k) w[row[k]] = nz[k];
    // Add diagonal shift
    w[c] += v;
    // Copy content of column to new matrix
    for (k = colind_new[c]; k < colind_new[c + 1]; ++k) nz[k] = w[row_new[k]];
    // Restore work vector
    for (k = colind[c]; k < colind[c + 1]; ++k) w[row[k]] = 0;
    w[c] = 0;
  }
}

} // namespace casadi
