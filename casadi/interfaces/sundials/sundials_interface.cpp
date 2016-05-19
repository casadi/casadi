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
        "Absolute tolerence  for the IVP solution"}},
      {"exact_jacobian",
       {OT_BOOL,
        "Use exact Jacobian information for the forward integration"}},
      {"exact_jacobianB",
       {OT_BOOL,
        "Use exact Jacobian information for the backward integration "
        "[default: equal to exact_jacobian]"}},
      {"upper_bandwidth",
       {OT_INT,
        "Upper band-width of banded Jacobian (estimations)"}},
      {"lower_bandwidth",
       {OT_INT,
        "Lower band-width of banded Jacobian (estimations)"}},
      {"linear_solver_type",
       {OT_STRING,
        "Type of iterative solver: USER_DEFINED|dense|banded|iterative"}},
      {"iterative_solver",
       {OT_STRING,
        "Iterative solver: GMRES|bcgstab|tfqmr"}},
      {"pretype",
       {OT_STRING,
        "Type of preconditioning: NONE|left|right|both"}},
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
        "Precondition an iterative solver"}},
      {"use_preconditionerB",
       {OT_BOOL,
        "Precondition an iterative solver for the backwards problem "
        "[default: equal to use_preconditioner]"}},
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
      {"finite_difference_fsens",
       {OT_BOOL,
        "Use finite differences to approximate the forward sensitivity equations "
        "(if AD is not available)"}},
      {"fsens_reltol",
       {OT_DOUBLE,
        "Relative tolerence for the forward sensitivity solution [default: equal to reltol]"}},
      {"fsens_abstol",
       {OT_DOUBLE,
        "Absolute tolerence for the forward sensitivity solution [default: equal to abstol]"}},
      {"fsens_scaling_factors",
       {OT_DOUBLEVECTOR,
        "Scaling factor for the components if finite differences is used"}},
      {"fsens_sensitiviy_parameters",
       {OT_INTVECTOR,
        "Specifies which components will be used when estimating the sensitivity equations"}},
      {"steps_per_checkpoint",
       {OT_INT,
        "Number of steps between two consecutive checkpoints"}},
      {"interpolation_type",
       {OT_STRING,
        "Type of interpolation for the adjoint sensitivities"}},
      {"upper_bandwidthB",
       {OT_INT,
        "Upper band-width of banded jacobians for backward integration "
        "[default: equal to upper_bandwidth]"}},
      {"lower_bandwidthB",
       {OT_INT,
        "lower band-width of banded jacobians for backward integration "
        "[default: equal to lower_bandwidth]"}},
      {"linear_solver_typeB",
       {OT_STRING,
        "Linear solver for backward integration"}},
      {"iterative_solverB",
       {OT_STRING,
        "Iterative solver for backward integration"}},
      {"pretypeB",
       {OT_STRING,
        "Preconditioner for backward integration"}},
      {"max_krylovB",
       {OT_INT,
        "Maximum krylov subspace size for backward integration"}},
      {"reltolB",
       {OT_DOUBLE,
        "Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]"}},
      {"abstolB",
       {OT_DOUBLE,
        "Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]"}},
      {"linear_solver",
       {OT_STRING,
        "A custom linear solver creator function [default: csparse]"}},
      {"linear_solver_options",
       {OT_DICT,
        "Options to be passed to the linear solver"}},
      {"linear_solverB",
       {OT_STRING,
        "A custom linear solver creator function for backwards integration "
        "[default: equal to linear_solver]"}},
      {"linear_solver_optionsB",
       {OT_DICT,
        "Options to be passed to the linear solver for backwards integration "
        "[default: equal to linear_solver_options]"}}
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
    exact_jac_ = true;
    max_num_steps_ = 10000;
    finite_difference_fsens_ = false;
    stop_at_end_ = true;
    use_precon_ = false;
    max_krylov_ = 10;
    linear_solver_ = "csparse";
    string linear_solver_type = "user_defined";
    string iterative_solver = "gmres";
    string pretype = "none";
    ubw_ = -1;
    lbw_ = -1;
    ubwB_ = -1;
    lbwB_ = -1;
    quad_err_con_ = false;
    interpolation_type_ = "hermite";
    steps_per_checkpoint_ = 20;
    disable_internal_warnings_ = false;
    max_multistep_order_ = 5;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="abstol") {
        abstol_ = op.second;
      } else if (op.first=="reltol") {
        reltol_ = op.second;
      } else if (op.first=="exact_jacobian") {
        exact_jac_ = op.second;
      } else if (op.first=="max_num_steps") {
        max_num_steps_ = op.second;
      } else if (op.first=="finite_difference_fsens") {
        finite_difference_fsens_ = op.second;
      } else if (op.first=="stop_at_end") {
        stop_at_end_ = op.second;
      } else if (op.first=="use_preconditioner") {
        use_precon_ = op.second;
      } else if (op.first=="max_krylov") {
        max_krylov_ = op.second;
      } else if (op.first=="linear_solver_type") {
        linear_solver_type = op.second.to_string();
      } else if (op.first=="iterative_solver") {
        iterative_solver = op.second.to_string();
      } else if (op.first=="pretype") {
        pretype = op.second.to_string();
      } else if (op.first=="linear_solver") {
        linear_solver_ = op.second.to_string();
      } else if (op.first=="linear_solver_options") {
        linear_solver_options_ = op.second;
      } else if (op.first=="upper_bandwidth") {
        ubw_ = op.second;
      } else if (op.first=="lower_bandwidth") {
        lbw_ = op.second;
      } else if (op.first=="upper_bandwidthB") {
        ubwB_ = op.second;
      } else if (op.first=="lower_bandwidthB") {
        lbwB_ = op.second;
      } else if (op.first=="quad_err_con") {
        quad_err_con_ = op.second;
      } else if (op.first=="interpolation_type") {
        interpolation_type_ = op.second.to_string();
      } else if (op.first=="steps_per_checkpoint") {
        steps_per_checkpoint_ = op.second;
      } else if (op.first=="disable_internal_warnings") {
        disable_internal_warnings_ = op.second;
      } else if (op.first=="max_multistep_order") {
        max_multistep_order_ = op.second;
      }
    }

    // Default dependent options
    exact_jacB_ = exact_jac_;
    fsens_abstol_ = abstol_;
    fsens_reltol_ = reltol_;
    abstolB_ = abstol_;
    reltolB_ = reltol_;
    use_preconB_ = use_precon_;
    max_krylovB_ = max_krylov_;
    std::string linear_solver_typeB = linear_solver_type;
    std::string iterative_solverB = iterative_solver;
    std::string pretypeB = pretype;
    linear_solverB_ = linear_solver_;
    linear_solver_optionsB_ = linear_solver_options_;

    // Read options again
    for (auto&& op : opts) {
      if (op.first=="exact_jacobianB") {
        exact_jacB_ = op.second;
      } else if (op.first=="fsens_abstol") {
        fsens_abstol_ = op.second;
      } else if (op.first=="fsens_reltol") {
        fsens_reltol_ = op.second;
      } else if (op.first=="abstolB") {
        abstolB_ = op.second;
      } else if (op.first=="reltolB") {
        reltolB_ = op.second;
      } else if (op.first=="use_preconditionerB") {
        use_preconB_ = op.second;
      } else if (op.first=="max_krylovB") {
        max_krylovB_ = op.second;
      } else if (op.first=="linear_solver_typeB") {
        linear_solver_typeB = op.second.to_string();
      } else if (op.first=="iterative_solverB") {
        iterative_solverB = op.second.to_string();
      } else if (op.first=="pretypeB") {
        pretypeB = op.second.to_string();
      } else if (op.first=="linear_solverB") {
        linear_solverB_ = op.second.to_string();
      } else if (op.first=="linear_solver_optionsB") {
        linear_solver_optionsB_ = op.second;
      }
    }

    // No Jacobian of g if g doesn't exist
    if (nrx_==0) {
      exact_jacB_ = false;
    }

    // Linear solver for forward integration
    if (linear_solver_type=="dense") {
      linsol_f_ = SD_DENSE;
    } else if (linear_solver_type=="banded") {
      linsol_f_ = SD_BANDED;
    } else if (linear_solver_type=="iterative") {
      linsol_f_ = SD_ITERATIVE;

      // Iterative solver
      if (iterative_solver=="gmres") {
        itsol_f_ = SD_GMRES;
      } else if (iterative_solver=="bcgstab") {
        itsol_f_ = SD_BCGSTAB;
      } else if (iterative_solver=="tfqmr") {
        itsol_f_ = SD_TFQMR;
      } else {
        casadi_error("Unknown iterative solver for forward integration: " + iterative_solver);
      }

      // Preconditioning type
      if (pretype=="none") {
        pretype_f_ = PREC_NONE;
      } else if (pretype=="left") {
        pretype_f_ = PREC_LEFT;
      } else if (pretype=="right") {
        pretype_f_ = PREC_RIGHT;
      } else if (pretype=="both") {
        pretype_f_ = PREC_BOTH;
      } else {
        casadi_error("Unknown preconditioning type for forward integration: " + pretype);
      }
    } else if (linear_solver_type=="user_defined") {
      linsol_f_ = SD_USER_DEFINED;
    } else {
      casadi_error("Unknown linear solver for forward integration: " + linear_solver_type);
    }

    // Linear solver for backward integration
    if (linear_solver_typeB=="dense") {
      linsol_g_ = SD_DENSE;
    } else if (linear_solver_typeB=="banded") {
      linsol_g_ = SD_BANDED;
    } else if (linear_solver_typeB=="iterative") {
      linsol_g_ = SD_ITERATIVE;

      // Iterative solver
      if (iterative_solverB=="gmres") {
        itsol_g_ = SD_GMRES;
      } else if (iterative_solverB=="bcgstab") {
        itsol_g_ = SD_BCGSTAB;
      } else if (iterative_solverB=="tfqmr") {
        itsol_g_ = SD_TFQMR;
      } else {
        casadi_error("Unknown sparse solver for backward integration: " + iterative_solverB);
      }

      // Preconditioning type
      if (pretypeB=="none") {
        pretype_g_ = PREC_NONE;
      } else if (pretypeB=="left") {
        pretype_g_ = PREC_LEFT;
      } else if (pretypeB=="right") {
        pretype_g_ = PREC_RIGHT;
      } else if (pretypeB=="both") {
        pretype_g_ = PREC_BOTH;
      } else {
        casadi_error("Unknown preconditioning type for backward integration: " + pretypeB);
      }
    } else if (linear_solver_typeB=="user_defined") {
      linsol_g_ = SD_USER_DEFINED;
    } else {
      casadi_error("Unknown linear solver for backward integration: " + iterative_solverB);
    }

    // Allocate work vectors
    alloc_w(np_, true); // p
    alloc_w(nrp_, true); // rp

    // Is exact Jacobian required?
    if (linsol_f_==SD_USER_DEFINED ||
      (linsol_f_==SD_ITERATIVE && use_precon_)) {
        casadi_assert_message(exact_jac_, "Exact Jacobian required");
    }

    // Is exact Jacobian required?
    if (nrx_>0) {
      if (linsol_g_==SD_USER_DEFINED ||
        (linsol_g_==SD_ITERATIVE && use_preconB_)) {
          casadi_assert_message(exact_jacB_, "Exact Jacobian required");
      }
    }
  }

  void SundialsInterface::init_jacF() {
    // Get sparsity
    const Sparsity& sp = get_function("jacF").sparsity_out(0);

    // Work vector
    alloc_w(sp.nnz(), true);

    // Create a linear solver
    set_function(linsol("linsolF", linear_solver_, sp, 1, linear_solver_options_));

    // Update bandwidths
    lbw_ = sp.bw_lower();
    ubw_ = sp.bw_upper();
  }

  void SundialsInterface::init_jacB() {
    // Get sparsity
    const Sparsity& sp = get_function("jacB").sparsity_out(0);

    // Work vector
    alloc_w(sp.nnz(), true);

    // Create a linear solver
    set_function(linsol("linsolB", linear_solverB_, sp, 1, linear_solver_optionsB_));

    // Update bandwidths
    lbwB_ = sp.bw_lower();
    ubwB_ = sp.bw_upper();
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

  void SundialsInterface::set_work(void* mem, const double**& arg, double**& res,
                                int*& iw, double*& w) const {
    auto m = static_cast<SundialsMemory*>(mem);

    // Set work in base classes
    Integrator::set_work(mem, arg, res, iw, w);

    // Work vectors
    m->p = w; w += np_;
    m->rp = w; w += nrp_;
    m->jac = w; w += get_function("jacF").nnz_out(0);
    if (nrx_>0) {
      m->jacB = w; w += get_function("jacB").nnz_out(0);
    }
  }

} // namespace casadi
