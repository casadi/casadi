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

using namespace std;
namespace casadi {

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
        "Coefficient in the nonlinear convergence test"}}
     }
  };

  void SundialsSimulator::init(const Dict& opts) {
    // Call the base class method
    Simulator::init(opts);

    // If sensitivity equations, make sure derivative_of_ is available
    casadi_assert(ns_==0 || !derivative_of_.is_null(),
      "Not implemented.");

    // Default options
    abstol_ = 1e-8;
    reltol_ = 1e-6;
    max_num_steps_ = 10000;
    stop_at_end_ = true;
    use_precon_ = true;
    max_krylov_ = 10;
    linear_solver_ = "qr";
    string newton_scheme = "direct";
    quad_err_con_ = false;
    string interpolation_type = "hermite";
    steps_per_checkpoint_ = 20;
    disable_internal_warnings_ = false;
    max_multistep_order_ = 5;
    second_order_correction_ = true;
    step0_ = 0;
    max_step_size_ = 0;
    max_order_ = 0;
    nonlin_conv_coeff_ = 0;

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

    // Get or create Jacobians and linear system solvers
    Function J;
    if (ns_==0) {
      J = getJ();
    } else {
      SundialsSimulator* d = derivative_of_.get<SundialsSimulator>();
      casadi_assert_dev(d!=nullptr);
      if (d->ns_==0) {
        J = d->get_function("jacF");
      } else {
        J = d->getJ();
      }
    }
    set_function(J, J.name(), true);
    alloc_w(J.nnz_out(0), true);

    // Allocate work vectors
    alloc_w(nu_, true); // u
    alloc_w(np_, true); // p
    alloc_w(2 * (nx_+nz_), true); // v1, v2

    // Allocate linear solvers
    linsolF_ = Linsol("linsolF", linear_solver_,
      get_function("jacF").sparsity_out(0), linear_solver_options_);
  }

  int SundialsSimulator::init_mem(void* mem) const {
    if (Simulator::init_mem(mem)) return 1;
    auto m = static_cast<SundialsSimMemory*>(mem);

    // Allocate n-vectors
    m->xz = N_VNew_Serial(nx_+nz_);

    m->mem_linsolF = linsolF_.checkout();

    return 0;
  }

  void SundialsSimulator::free_mem(void *mem) const {
    Simulator::free_mem(mem);
    auto m = static_cast<SundialsSimMemory*>(mem);

    linsolF_.release(m->mem_linsolF);
  }

  void SundialsSimulator::reset(SimulatorMemory* mem, double t, const double* x, const double* u,
      const double* z, const double* p, double* y) const {
    auto m = static_cast<SundialsSimMemory*>(mem);
    // Update time
    m->t = t;
    // Set parameters
    casadi_copy(p, np_, m->p);
    // Set controls
    casadi_copy(u, nu_, m->u);
    // Set the state
    casadi_copy(x, nx_, NV_DATA_S(m->xz));
    casadi_copy(z, nz_, NV_DATA_S(m->xz) + nx_);
  }

  SundialsSimMemory::SundialsSimMemory() {
    this->xz  = nullptr;
  }

  SundialsSimMemory::~SundialsSimMemory() {
    if (this->xz) N_VDestroy_Serial(this->xz);
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
    print("FORWARD INTEGRATION:\n");
    print("Number of steps taken by SUNDIALS: %ld\n", m->nsteps);
    print("Number of calls to the user’s f function: %ld\n", m->nfevals);
    print("Number of calls made to the linear solver setup function: %ld\n", m->nlinsetups);
    print("Number of error test failures: %ld\n", m->netfails);
    print("Method order used on the last internal step: %d\n", m->qlast);
    print("Method order to be used on the next internal step: %d\n", m->qcur);
    print("Actual value of initial step size: %g\n", m->hinused);
    print("Step size taken on the last internal step: %g\n", m->hlast);
    print("Step size to be attempted on the next internal step: %g\n", m->hcur);
    print("Current internal time reached: %g\n");
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
    m->u = w; w += nu_;
    m->p = w; w += np_;
    m->v1 = w; w += nx_ + nz_;
    m->v2 = w; w += nx_ + nz_;
    m->jac = w; w += get_function("jacF").nnz_out(0);
  }

} // namespace casadi
