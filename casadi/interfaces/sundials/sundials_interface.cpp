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

  SundialsInterface::SundialsInterface(const std::string& name, const XProblem& dae)
    : Integrator(name, dae) {
    addOption("max_num_steps",               OT_INT,          10000,
              "Maximum number of integrator steps");
    addOption("reltol",                      OT_DOUBLE,             1e-6,
              "Relative tolerence for the IVP solution");
    addOption("abstol",                      OT_DOUBLE,             1e-8,
              "Absolute tolerence  for the IVP solution");
    addOption("exact_jacobian",              OT_BOOL,          true,
              "Use exact Jacobian information for the forward integration");
    addOption("exact_jacobianB",             OT_BOOL,          GenericType(),
              "Use exact Jacobian information for the backward integration "
              "[default: equal to exact_jacobian]");
    addOption("upper_bandwidth",             OT_INT,          GenericType(),
              "Upper band-width of banded Jacobian (estimations)");
    addOption("lower_bandwidth",             OT_INT,          GenericType(),
              "Lower band-width of banded Jacobian (estimations)");
    addOption("linear_solver_type",          OT_STRING,           "dense", "");
    addOption("iterative_solver",            OT_STRING,           "gmres", "");
    addOption("pretype",                     OT_STRING,           "none", "");
    addOption("max_krylov",                  OT_INT,          10,
              "Maximum Krylov subspace size");
    addOption("sensitivity_method",          OT_STRING,           "simultaneous", "");
    addOption("max_multistep_order",         OT_INT,          5);
    addOption("use_preconditioner",          OT_BOOL,          false,
              "Precondition an iterative solver");
    addOption("use_preconditionerB",         OT_BOOL,          GenericType(),
              "Precondition an iterative solver for the backwards problem "
              "[default: equal to use_preconditioner]");
    addOption("stop_at_end",                 OT_BOOL,          true,
              "Stop the integrator at the end of the interval");
    addOption("disable_internal_warnings",        OT_BOOL,             false,
              "Disable SUNDIALS internal warning messages");

    // Quadratures
    addOption("quad_err_con",                OT_BOOL,          false,
              "Should the quadratures affect the step size control");

    // Forward sensitivity problem
    addOption("fsens_err_con",               OT_BOOL,          true,
              "include the forward sensitivities in all error controls");
    addOption("finite_difference_fsens",     OT_BOOL,          false,
              "Use finite differences to approximate the forward sensitivity equations "
              "(if AD is not available)");
    addOption("fsens_reltol",                OT_DOUBLE,             GenericType(),
              "Relative tolerence for the forward sensitivity solution [default: equal to reltol]");
    addOption("fsens_abstol",                OT_DOUBLE,             GenericType(),
              "Absolute tolerence for the forward sensitivity solution [default: equal to abstol]");
    addOption("fsens_scaling_factors",       OT_DOUBLEVECTOR,       GenericType(),
              "Scaling factor for the components if finite differences is used");
    addOption("fsens_sensitiviy_parameters", OT_INTVECTOR,    GenericType(),
              "Specifies which components will be used when estimating the sensitivity equations");

    // Adjoint sensivity problem
    addOption("steps_per_checkpoint",        OT_INT,          20,
              "Number of steps between two consecutive checkpoints");
    addOption("interpolation_type",          OT_STRING,           "hermite",
              "Type of interpolation for the adjoint sensitivities");
    addOption("upper_bandwidthB",            OT_INT,          GenericType(),
              "Upper band-width of banded jacobians for backward integration "
              "[default: equal to upper_bandwidth]");
    addOption("lower_bandwidthB",            OT_INT,          GenericType(),
              "lower band-width of banded jacobians for backward integration "
              "[default: equal to lower_bandwidth]");
    addOption("linear_solver_typeB",         OT_STRING,           GenericType(), "");
    addOption("iterative_solverB",           OT_STRING,           GenericType(), "");
    addOption("pretypeB",                    OT_STRING,           GenericType(), "");
    addOption("max_krylovB",                 OT_INT,          GenericType(),
              "Maximum krylov subspace size");
    addOption("reltolB",                     OT_DOUBLE,             GenericType(),
              "Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]");
    addOption("abstolB",                     OT_DOUBLE,             GenericType(),
              "Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]");
    addOption("linear_solver",               OT_STRING,     GenericType(),
              "A custom linear solver creator function");
    addOption("linear_solver_options",       OT_DICT,       GenericType(),
              "Options to be passed to the linear solver");
    addOption("linear_solverB",              OT_STRING,     GenericType(),
              "A custom linear solver creator function for backwards integration "
              "[default: equal to linear_solver]");
    addOption("linear_solver_optionsB",      OT_DICT,       GenericType(),
              "Options to be passed to the linear solver for backwards integration "
              "[default: equal to linear_solver_options]");
  }

  SundialsInterface::~SundialsInterface() {
  }

  void SundialsInterface::init(const Dict& opts) {
    // Call the base class method
    Integrator::init(opts);

    // Default options
    abstol_ = 1e-8;
    reltol_ = 1e-6;
    exact_jacobian_ = true;
    max_num_steps_ = 10000;
    finite_difference_fsens_ = false;
    stop_at_end_ = true;
    use_preconditioner_ = false;
    max_krylov_ = 10;
    string linear_solver_type = "dense";
    string iterative_solver = "gmres";
    string pretype = "none";
    string linear_solver;
    Dict linear_solver_options;
    upper_bandwidth_ = -1;
    lower_bandwidth_ = -1;
    upper_bandwidthB_ = -1;
    lower_bandwidthB_ = -1;
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
        exact_jacobian_ = op.second;
      } else if (op.first=="max_num_steps") {
        max_num_steps_ = op.second;
      } else if (op.first=="finite_difference_fsens") {
        finite_difference_fsens_ = op.second;
      } else if (op.first=="stop_at_end") {
        stop_at_end_ = op.second;
      } else if (op.first=="use_preconditioner") {
        use_preconditioner_ = op.second;
      } else if (op.first=="max_krylov") {
        max_krylov_ = op.second;
      } else if (op.first=="linear_solver_type") {
        linear_solver_type = op.second.to_string();
      } else if (op.first=="iterative_solver") {
        iterative_solver = op.second.to_string();
      } else if (op.first=="pretype") {
        pretype = op.second.to_string();
      } else if (op.first=="linear_solver") {
        linear_solver = op.second.to_string();
      } else if (op.first=="linear_solver_options") {
        linear_solver_options = op.second;
      } else if (op.first=="upper_bandwidth") {
        upper_bandwidth_ = op.second;
      } else if (op.first=="lower_bandwidth") {
        lower_bandwidth_ = op.second;
      } else if (op.first=="upper_bandwidthB") {
        upper_bandwidthB_ = op.second;
      } else if (op.first=="lower_bandwidthB") {
        lower_bandwidthB_ = op.second;
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
    exact_jacobianB_ = exact_jacobian_;
    fsens_abstol_ = abstol_;
    fsens_reltol_ = reltol_;
    abstolB_ = abstol_;
    reltolB_ = reltol_;
    use_preconditionerB_ = use_preconditioner_;
    max_krylovB_ = max_krylov_;
    std::string linear_solver_typeB = linear_solver_type;
    std::string iterative_solverB = iterative_solver;
    std::string pretypeB = pretype;
    string linear_solverB = linear_solver;
    Dict linear_solver_optionsB = linear_solver_options;

    // Read options again
    for (auto&& op : opts) {
      if (op.first=="exact_jacobianB") {
        exact_jacobianB_ = op.second;
      } else if (op.first=="fsens_abstol") {
        fsens_abstol_ = op.second;
      } else if (op.first=="fsens_reltol") {
        fsens_reltol_ = op.second;
      } else if (op.first=="abstolB") {
        abstolB_ = op.second;
      } else if (op.first=="reltolB") {
        reltolB_ = op.second;
      } else if (op.first=="use_preconditionerB") {
        use_preconditionerB_ = op.second;
      } else if (op.first=="max_krylovB") {
        max_krylovB_ = op.second;
      } else if (op.first=="linear_solver_typeB") {
        linear_solver_typeB = op.second.to_string();
      } else if (op.first=="iterative_solverB") {
        iterative_solverB = op.second.to_string();
      } else if (op.first=="pretypeB") {
        pretypeB = op.second.to_string();
      } else if (op.first=="linear_solverB") {
        linear_solverB = op.second.to_string();
      } else if (op.first=="linear_solver_optionsB") {
        linear_solver_optionsB = op.second;
      }
    }

    // No Jacobian of g if g doesn't exist
    if (g_.isNull()) {
      exact_jacobianB_ = false;
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

    // Create a Jacobian if requested
    if (exact_jacobian_) {
      jac_ = getJac();
      alloc(jac_);
      alloc_w(jac_.nnz_out(0), true);
    }

    if (!jac_.isNull()) {
      casadi_assert_message(jac_.size2_out(0)==jac_.size1_out(0),
                            "SundialsInterface::init: the jacobian of the forward problem must "
                            "be square but got " << jac_.sparsity_out(0).dim());

      casadi_assert_message(!jac_.sparsity_out(0).is_singular(),
                            "SundialsInterface::init: singularity - the jacobian of the forward "
                            "problem is structurally rank-deficient. sprank(J)="
                            << sprank(jac_.sparsity_out(0)) << " (in stead of "<< jac_.size2_out(0)
                            << ")");
    }

    // Create a backwards Jacobian if requested
    if (exact_jacobianB_ && !g_.isNull()) jacB_ = getJacB();

    if (!jacB_.isNull()) {
      alloc(jacB_);
      alloc_w(jacB_.nnz_out(0), true);
      casadi_assert_message(jacB_.size2_out(0)==jacB_.size1_out(0),
                            "SundialsInterface::init: the jacobian of the backward problem must be "
                            "square but got " << jacB_.sparsity_out(0).dim());

      casadi_assert_message(!jacB_.sparsity_out(0).is_singular(),
                            "SundialsInterface::init: singularity - the jacobian of the backward"
                            " problem is structurally rank-deficient. sprank(J)="
                            << sprank(jacB_.sparsity_out(0)) << " (instead of "
                            << jacB_.size2_out(0) << ")");
    }

    // Create a linear solver
    if (!linear_solver.empty() && !jac_.isNull()) {
      linsol_ = linsol("linsol", linear_solver, jac_.sparsity_out(0),
                       1, linear_solver_options);
      alloc(linsol_);
    }

    // Create a linear solver
    if (!linear_solverB.empty() && !jacB_.isNull()) {
      linsolB_ = linsol("linsolB", linear_solverB, jacB_.sparsity_out(0),
                        1, linear_solver_optionsB);
      alloc(linsolB_);
    }

    // Allocate temporary memory
    //alloc_w(np_, true); // p_
    //alloc_w(nrp_, true); // rp_
  }

  void SundialsInterface::init_memory(Memory& mem) const {
    Integrator::init_memory(mem);
    SundialsMemory& m = dynamic_cast<SundialsMemory&>(mem);

    // Allocate n-vectors
    m.xz = N_VNew_Serial(nx_+nz_);
    m.q = N_VNew_Serial(nq_);
    m.rxz = N_VNew_Serial(nrx_+nrz_);
    m.rq = N_VNew_Serial(nrq_);

    // Allocate memory
    m.p.resize(np_);
    m.rp.resize(nrp_);
  }

  void SundialsInterface::reset(IntegratorMemory& mem, double t, const double* x,
                                const double* z, const double* p) const {
    SundialsMemory& m = dynamic_cast<SundialsMemory&>(mem);

    // Update time
    m.t = t;

    // Set parameters
    casadi_copy(p, np_, get_ptr(m.p));

    // Set the state
    casadi_copy(x, nx_, NV_DATA_S(m.xz));
    casadi_copy(z, nz_, NV_DATA_S(m.xz)+nx_);

    // Reset summation states
    N_VConst(0., m.q);
  }

  void SundialsInterface::resetB(IntegratorMemory& mem, double t, const double* rx,
                                 const double* rz, const double* rp) const {
    SundialsMemory& m = dynamic_cast<SundialsMemory&>(mem);

    // Update time
    m.t = t;

    // Set parameters
    casadi_copy(rp, nrp_, get_ptr(m.rp));

    // Get the backward state
    casadi_copy(rx, nrx_, NV_DATA_S(m.rxz));

    // Reset summation states
    N_VConst(0., m.rq);
  }

  std::pair<int, int> SundialsInterface::getBandwidth() const {
    std::pair<int, int> bw;

    // Get upper bandwidth
    if (upper_bandwidth_>=0) {
      bw.first = upper_bandwidth_;
    } else {
      casadi_assert_message(!jac_.isNull(),
                            "\"upper_bandwidth\" has not been set and cannot be "
                            "detected since exact Jacobian is not available.");
      bw.first = jac_.sparsity_out(0).bw_upper();
    }

    // Get lower bandwidth
    if (lower_bandwidth_>=0) {
      bw.second = lower_bandwidth_;
    } else {
      casadi_assert_message(!jac_.isNull(),
                            "\"lower_bandwidth\" has not been set and cannot be "
                            "detected since exact Jacobian is not available.");
      bw.second = jac_.sparsity_out(0).bw_lower();
    }

    return bw;
  }

  std::pair<int, int> SundialsInterface::getBandwidthB() const {
    std::pair<int, int> bw;

    // Get upper bandwidth
    if (upper_bandwidthB_>=0) {
      bw.first = upper_bandwidthB_;
    } else {
      casadi_assert_message(!jacB_.isNull(),
                            "\"upper_bandwidthB\" has not been set and cannot be "
                            "detected since exact Jacobian for backward problem "
                            "is not available.");
      bw.first = jacB_.sparsity_out(0).bw_upper();
    }

    // Get lower bandwidth
    if (lower_bandwidthB_>=0) {
      bw.second = lower_bandwidthB_;
    } else {
      casadi_assert_message(!jacB_.isNull(),
                            "\"lower_bandwidthB\" has not been set and cannot be "
                            "detected since exact Jacobian for backward problem "
                            "is not available.");
      bw.second = jacB_.sparsity_out(0).bw_lower();
    }

    return bw;
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

  Dict SundialsMemory::get_stats() const {
    Dict stats = IntegratorMemory::get_stats();
    // Counters, forward problem
    stats["nsteps"] = static_cast<int>(nsteps);
    stats["nfevals"] = static_cast<int>(nfevals);
    stats["nlinsetups"] = static_cast<int>(nlinsetups);
    stats["netfails"] = static_cast<int>(netfails);
    stats["qlast"] = qlast;
    stats["qcur"] = qcur;
    stats["hinused"] = hinused;
    stats["hlast"] = hlast;
    stats["hcur"] = hcur;
    stats["tcur"] = tcur;

    // Counters, backward problem
    stats["nstepsB"] = static_cast<int>(nstepsB);
    stats["nfevalsB"] = static_cast<int>(nfevalsB);
    stats["nlinsetupsB"] = static_cast<int>(nlinsetupsB);
    stats["netfailsB"] = static_cast<int>(netfailsB);
    stats["qlastB"] = qlastB;
    stats["qcurB"] = qcurB;
    stats["hinusedB"] = hinusedB;
    stats["hlastB"] = hlastB;
    stats["hcurB"] = hcurB;
    stats["tcurB"] = tcurB;

    // Timers
    stats["t_res"] = t_res;
    stats["t_fres"] = t_fres;
    stats["t_jac"] = t_jac;
    stats["t_jacB"] = t_jacB;
    stats["t_lsolve"] = t_lsolve;
    stats["t_lsetup_jac"] = t_lsetup_jac;
    stats["t_lsetup_fac"] = t_lsetup_fac;
    return stats;
  }

} // namespace casadi


