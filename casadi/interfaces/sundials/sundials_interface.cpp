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


#include "sundials_interface.hpp"

#include "casadi/core/std_vector_tools.hpp"

INPUTSCHEME(IntegratorInput)
OUTPUTSCHEME(IntegratorOutput)

using namespace std;
namespace casadi {

  SundialsInterface::SundialsInterface(const std::string& name, const Function& dae)
    : Integrator(name, dae) {
  }

  SundialsInterface::~SundialsInterface() {
  }

  Options SundialsInterface::options_
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
      {"use_iterative_solver",
       {OT_BOOL,
        "Use iterative solver as opposed to a direct solver"}},
      {"iterative_solver",
       {OT_STRING,
        "Iterative solver: GMRES|bcgstab|tfqmr"}},
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
        "Stop the integrator at the end of the interval"}},
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
        "A custom linear solver creator function [default: csparse]"}},
      {"linear_solver_options",
       {OT_DICT,
        "Options to be passed to the linear solver"}}
     }
  };

  void SundialsInterface::init(const Dict& opts) {
    // Call the base class method
    Integrator::init(opts);

    // If sensitivity equations, make sure derivative_of_ is available
    casadi_assert_message(ns_==0 || !derivative_of_.is_null(),
      "Not implemented.");

    // Default options
    abstol_ = 1e-8;
    reltol_ = 1e-6;
    max_num_steps_ = 10000;
    stop_at_end_ = true;
    use_precon_ = true;
    max_krylov_ = 10;
    linear_solver_ = "csparse";
    use_iterative_solver_ = false;
    string iterative_solver = "gmres";
    quad_err_con_ = false;
    string interpolation_type = "hermite";
    steps_per_checkpoint_ = 20;
    disable_internal_warnings_ = false;
    max_multistep_order_ = 5;

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
      } else if (op.first=="use_iterative_solver") {
        use_iterative_solver_ = op.second;
      } else if (op.first=="iterative_solver") {
        iterative_solver = op.second.to_string();
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
      }
    }

    // Linear solver for forward integration
    if (use_iterative_solver_) {
      if (iterative_solver=="gmres") {
        itsol_ = SD_GMRES;
      } else if (iterative_solver=="bcgstab") {
        itsol_ = SD_BCGSTAB;
      } else if (iterative_solver=="tfqmr") {
        itsol_ = SD_TFQMR;
      } else {
        casadi_error("Unknown iterative solver: " + iterative_solver);
      }
    }

    // interpolation_type
    if (interpolation_type=="hermite") {
      interp_ = SD_HERMITE;
    } else if (interpolation_type=="polynomial") {
      interp_ = SD_POLYNOMIAL;
    } else {
      casadi_error("Unknown interpolation type: " + interpolation_type);
    }

    // Get or create Jacobians and linear system solvers
    for (bool backward : {false, true}) {
      // Skip backward?
      if (backward && nrx_==0) continue;

      // Get Jacobian function
      Function J, solver;
      if (ns_==0) {
        J = getJ(backward);
      } else {
        SundialsInterface* d = derivative_of_.get<SundialsInterface>();
        casadi_assert(d!=0);
        if (d->ns_==0) {
          J = d->get_function(backward ? "jacB" : "jacF");
          solver = d->get_function(backward ? "linsolB" : "linsolF");
        } else {
          J = d->getJ(backward);
        }
        // Augmented Jacobian-times-vector function
        Function jtimes = d->get_function(backward ? "jtimesB" : "jtimesF");
        set_function(jtimes.forward(1), backward ? "jacBv" : "jacFv");
      }
      set_function(J);
      alloc_w(J.nnz_out(0), true);

      // Allocate a linear solver
      if (solver.is_null()) {
        solver = linsol(backward ? "linsolB" : "linsolF", linear_solver_,
                        J.sparsity_out(0), 0, linear_solver_options_);
      }
      set_function(solver);
    }

    // Allocate work vectors
    alloc_w(np_, true); // p
    alloc_w(nrp_, true); // rp
    if (ns_>0) {
      alloc_w(max(nz_, nrz_), true); // ztmp
    }
  }

  void SundialsInterface::init_memory(void* mem) const {
    Integrator::init_memory(mem);
    auto m = static_cast<SundialsMemory*>(mem);

    // Allocate n-vectors
    m->xz = N_VNew_Serial(nx_+nz_);
    m->q = N_VNew_Serial(nq_);
    m->rxz = N_VNew_Serial(nrx_+nrz_);
    m->rq = N_VNew_Serial(nrq_);
  }

  void SundialsInterface::reset(IntegratorMemory* mem, double t, const double* x,
                                const double* z, const double* p) const {
    auto m = static_cast<SundialsMemory*>(mem);

    // Update time
    m->t = t;

    // Set parameters
    casadi_copy(p, np_, m->p);

    // Set the state
    casadi_copy(x, nx_, NV_DATA_S(m->xz));
    casadi_copy(z, nz_, NV_DATA_S(m->xz)+nx_);

    // Reset summation states
    N_VConst(0., m->q);
  }

  void SundialsInterface::resetB(IntegratorMemory* mem, double t, const double* rx,
                                 const double* rz, const double* rp) const {
    auto m = static_cast<SundialsMemory*>(mem);

    // Update time
    m->t = t;

    // Set parameters
    casadi_copy(rp, nrp_, m->rp);

    // Get the backward state
    casadi_copy(rx, nrx_, NV_DATA_S(m->rxz));

    // Reset summation states
    N_VConst(0., m->rq);
  }

  SundialsMemory::SundialsMemory() {
    this->xz  = 0;
    this->q = 0;
    this->rxz = 0;
    this->rq = 0;
    this->first_callB = true;
  }

  SundialsMemory::~SundialsMemory() {
    if (this->xz) N_VDestroy_Serial(this->xz);
    if (this->q) N_VDestroy_Serial(this->q);
    if (this->rxz) N_VDestroy_Serial(this->rxz);
    if (this->rq) N_VDestroy_Serial(this->rq);
  }

  Dict SundialsInterface::get_stats(void* mem) const {
    Dict stats = Integrator::get_stats(mem);
    auto m = static_cast<SundialsMemory*>(mem);

    // Counters, forward problem
    stats["nsteps"] = static_cast<int>(m->nsteps);
    stats["nfevals"] = static_cast<int>(m->nfevals);
    stats["nlinsetups"] = static_cast<int>(m->nlinsetups);
    stats["netfails"] = static_cast<int>(m->netfails);
    stats["qlast"] = m->qlast;
    stats["qcur"] = m->qcur;
    stats["hinused"] = m->hinused;
    stats["hlast"] = m->hlast;
    stats["hcur"] = m->hcur;
    stats["tcur"] = m->tcur;

    // Counters, backward problem
    stats["nstepsB"] = static_cast<int>(m->nstepsB);
    stats["nfevalsB"] = static_cast<int>(m->nfevalsB);
    stats["nlinsetupsB"] = static_cast<int>(m->nlinsetupsB);
    stats["netfailsB"] = static_cast<int>(m->netfailsB);
    stats["qlastB"] = m->qlastB;
    stats["qcurB"] = m->qcurB;
    stats["hinusedB"] = m->hinusedB;
    stats["hlastB"] = m->hlastB;
    stats["hcurB"] = m->hcurB;
    stats["tcurB"] = m->tcurB;
    return stats;
  }

  void SundialsInterface::print_stats(IntegratorMemory* mem, ostream &stream) const {
    auto m = to_mem(mem);
    stream << "FORWARD INTEGRATION:" << endl;
    stream << "Number of steps taken by SUNDIALS: " << m->nsteps << endl;
    stream << "Number of calls to the user’s f function: " << m->nfevals << endl;
    stream << "Number of calls made to the linear solver setup function: "
           << m->nlinsetups << endl;
    stream << "Number of error test failures: " << m->netfails << endl;
    stream << "Method order used on the last internal step: "  << m->qlast << endl;
    stream << "Method order to be used on the next internal step: " << m->qcur << endl;
    stream << "Actual value of initial step size: " << m->hinused << endl;
    stream << "Step size taken on the last internal step: " << m->hlast << endl;
    stream << "Step size to be attempted on the next internal step: " << m->hcur << endl;
    stream << "Current internal time reached: " << m->tcur << endl;
    stream << "Number of checkpoints stored: " << m->ncheck << endl;
    if (nrx_>0) {
      stream << "BACKWARD INTEGRATION:" << endl;
      stream << "Number of steps taken by SUNDIALS: " << m->nstepsB << endl;
      stream << "Number of calls to the user’s f function: " << m->nfevalsB << endl;
      stream << "Number of calls made to the linear solver setup function: "
             << m->nlinsetupsB << endl;
      stream << "Number of error test failures: " << m->netfailsB << endl;
      stream << "Method order used on the last internal step: "  << m->qlastB << endl;
      stream << "Method order to be used on the next internal step: " << m->qcurB << endl;
      stream << "Actual value of initial step size: " << m->hinusedB << endl;
      stream << "Step size taken on the last internal step: " << m->hlastB << endl;
      stream << "Step size to be attempted on the next internal step: " << m->hcurB << endl;
      stream << "Current internal time reached: " << m->tcurB << endl;
    }
    stream << endl;
  }

  void SundialsInterface::set_work(void* mem, const double**& arg, double**& res,
                                int*& iw, double*& w) const {
    auto m = static_cast<SundialsMemory*>(mem);

    // Set work in base classes
    Integrator::set_work(mem, arg, res, iw, w);

    // Work vectors
    m->p = w; w += np_;
    m->rp = w; w += nrp_;
    if (ns_>0) {
      m->ztmp = w; w += max(nz_, nrz_);
    }
    m->jac = w; w += get_function("jacF").nnz_out(0);
    if (nrx_>0) {
      m->jacB = w; w += get_function("jacB").nnz_out(0);
    }
  }

} // namespace casadi
