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

INPUTSCHEME(IvpsolInput)
OUTPUTSCHEME(IvpsolOutput)

using namespace std;
namespace casadi {

  SundialsInterface::SundialsInterface(const std::string& name, const XProblem& dae)
    : Ivpsol(name, dae) {
    addOption("max_num_steps",               OT_INTEGER,          10000,
              "Maximum number of integrator steps");
    addOption("reltol",                      OT_REAL,             1e-6,
              "Relative tolerence for the IVP solution");
    addOption("abstol",                      OT_REAL,             1e-8,
              "Absolute tolerence  for the IVP solution");
    addOption("exact_jacobian",              OT_BOOLEAN,          true,
              "Use exact Jacobian information for the forward integration");
    addOption("exact_jacobianB",             OT_BOOLEAN,          GenericType(),
              "Use exact Jacobian information for the backward integration "
              "[default: equal to exact_jacobian]");
    addOption("upper_bandwidth",             OT_INTEGER,          GenericType(),
              "Upper band-width of banded Jacobian (estimations)");
    addOption("lower_bandwidth",             OT_INTEGER,          GenericType(),
              "Lower band-width of banded Jacobian (estimations)");
    addOption("linear_solver_type",          OT_STRING,           "dense",
              "", "user_defined|dense|banded|iterative");
    addOption("iterative_solver",            OT_STRING,           "gmres",
              "", "gmres|bcgstab|tfqmr");
    addOption("pretype",                     OT_STRING,           "none",
              "", "none|left|right|both");
    addOption("max_krylov",                  OT_INTEGER,          10,
              "Maximum Krylov subspace size");
    addOption("sensitivity_method",          OT_STRING,           "simultaneous", "",
              "simultaneous|staggered");
    addOption("max_multistep_order",         OT_INTEGER,          5);
    addOption("use_preconditioner",          OT_BOOLEAN,          false,
              "Precondition an iterative solver");
    addOption("use_preconditionerB",         OT_BOOLEAN,          GenericType(),
              "Precondition an iterative solver for the backwards problem "
              "[default: equal to use_preconditioner]");
    addOption("stop_at_end",                 OT_BOOLEAN,          true,
              "Stop the integrator at the end of the interval");

    // Quadratures
    addOption("quad_err_con",                OT_BOOLEAN,          false,
              "Should the quadratures affect the step size control");

    // Forward sensitivity problem
    addOption("fsens_err_con",               OT_BOOLEAN,          true,
              "include the forward sensitivities in all error controls");
    addOption("finite_difference_fsens",     OT_BOOLEAN,          false,
              "Use finite differences to approximate the forward sensitivity equations "
              "(if AD is not available)");
    addOption("fsens_reltol",                OT_REAL,             GenericType(),
              "Relative tolerence for the forward sensitivity solution [default: equal to reltol]");
    addOption("fsens_abstol",                OT_REAL,             GenericType(),
              "Absolute tolerence for the forward sensitivity solution [default: equal to abstol]");
    addOption("fsens_scaling_factors",       OT_REALVECTOR,       GenericType(),
              "Scaling factor for the components if finite differences is used");
    addOption("fsens_sensitiviy_parameters", OT_INTEGERVECTOR,    GenericType(),
              "Specifies which components will be used when estimating the sensitivity equations");

    // Adjoint sensivity problem
    addOption("steps_per_checkpoint",        OT_INTEGER,          20,
              "Number of steps between two consecutive checkpoints");
    addOption("interpolation_type",          OT_STRING,           "hermite",
              "Type of interpolation for the adjoint sensitivities", "hermite|polynomial");
    addOption("upper_bandwidthB",            OT_INTEGER,          GenericType(),
              "Upper band-width of banded jacobians for backward integration "
              "[default: equal to upper_bandwidth]");
    addOption("lower_bandwidthB",            OT_INTEGER,          GenericType(),
              "lower band-width of banded jacobians for backward integration "
              "[default: equal to lower_bandwidth]");
    addOption("linear_solver_typeB",         OT_STRING,           GenericType(),
              "", "user_defined|dense|banded|iterative");
    addOption("iterative_solverB",           OT_STRING,           GenericType(),
              "", "gmres|bcgstab|tfqmr");
    addOption("pretypeB",                    OT_STRING,           GenericType(),
              "", "none|left|right|both");
    addOption("max_krylovB",                 OT_INTEGER,          GenericType(),
              "Maximum krylov subspace size");
    addOption("reltolB",                     OT_REAL,             GenericType(),
              "Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]");
    addOption("abstolB",                     OT_REAL,             GenericType(),
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
    xz_  = 0;
    q_ = 0;
    rxz_ = 0;
    rq_ = 0;
  }

  SundialsInterface::~SundialsInterface() {
    if (xz_) { N_VDestroy_Serial(xz_); xz_ = 0; }
    if (q_) { N_VDestroy_Serial(q_); q_ = 0; }
    if (rxz_) { N_VDestroy_Serial(rxz_); rxz_ = 0; }
    if (rq_) { N_VDestroy_Serial(rq_); rq_ = 0; }
  }

  void SundialsInterface::init() {
    // Call the base class method
    Ivpsol::init();

    // Allocate n-vectors
    xz_ = N_VNew_Serial(nx_+nz_);
    q_ = N_VNew_Serial(nq_);
    rxz_ = N_VNew_Serial(nrx_+nrz_);
    rq_ = N_VNew_Serial(nrq_);

    // Read options
    abstol_ = option("abstol");
    reltol_ = option("reltol");
    exact_jacobian_ = option("exact_jacobian");
    exact_jacobianB_ = hasSetOption("exact_jacobianB") ?
      option("exact_jacobianB") && !g_.isNull() : exact_jacobian_;
    max_num_steps_ = option("max_num_steps");
    finite_difference_fsens_ = option("finite_difference_fsens");
    fsens_abstol_ =
      hasSetOption("fsens_abstol") ? static_cast<double>(option("fsens_abstol")) : abstol_;
    fsens_reltol_ =
      hasSetOption("fsens_reltol") ? static_cast<double>(option("fsens_reltol")) : reltol_;
    abstolB_ = hasSetOption("abstolB") ? static_cast<double>(option("abstolB")) : abstol_;
    reltolB_ = hasSetOption("reltolB") ? static_cast<double>(option("reltolB")) : reltol_;
    stop_at_end_ = option("stop_at_end");
    use_preconditioner_ = option("use_preconditioner");
    use_preconditionerB_ =  hasSetOption("use_preconditionerB") ?
      static_cast<bool>(option("use_preconditionerB")): use_preconditioner_;
    max_krylov_ = option("max_krylov");
    max_krylovB_ =
      hasSetOption("max_krylovB") ? static_cast<int>(option("max_krylovB")): max_krylov_;

    // Linear solver for forward integration
    if (option("linear_solver_type")=="dense") {
      linsol_f_ = SD_DENSE;
    } else if (option("linear_solver_type")=="banded") {
      linsol_f_ = SD_BANDED;
    } else if (option("linear_solver_type")=="iterative") {
      linsol_f_ = SD_ITERATIVE;

      // Iterative solver
      if (option("iterative_solver")=="gmres") {
        itsol_f_ = SD_GMRES;
      } else if (option("iterative_solver")=="bcgstab") {
        itsol_f_ = SD_BCGSTAB;
      } else if (option("iterative_solver")=="tfqmr") {
        itsol_f_ = SD_TFQMR;
      } else {
        throw CasadiException("Unknown sparse solver for forward integration");
      }

      // Preconditioning type
      if (option("pretype")=="none")               pretype_f_ = PREC_NONE;
      else if (option("pretype")=="left")          pretype_f_ = PREC_LEFT;
      else if (option("pretype")=="right")         pretype_f_ = PREC_RIGHT;
      else if (option("pretype")=="both")          pretype_f_ = PREC_BOTH;
      else
        throw CasadiException("Unknown preconditioning type for forward integration");
    } else if (option("linear_solver_type")=="user_defined") {
      linsol_f_ = SD_USER_DEFINED;
    } else {
      throw CasadiException("Unknown linear solver for forward integration");
    }


    std::string linear_solver_typeB = hasSetOption("linear_solver_typeB") ?
      option("linear_solver_typeB") : option("linear_solver_type");
    std::string iterative_solverB = hasSetOption("iterative_solverB") ?
      option("iterative_solverB") : option("iterative_solver");
    std::string pretypeB = hasSetOption("pretypeB") ? option("pretypeB"): option("pretype");

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
        throw CasadiException("Unknown sparse solver for backward integration");
      }

      // Preconditioning type
      if (pretypeB=="none")               pretype_g_ = PREC_NONE;
      else if (pretypeB=="left")          pretype_g_ = PREC_LEFT;
      else if (pretypeB=="right")         pretype_g_ = PREC_RIGHT;
      else if (pretypeB=="both")          pretype_g_ = PREC_BOTH;
      else
        throw CasadiException("Unknown preconditioning type for backward integration");
    } else if (linear_solver_typeB=="user_defined") {
      linsol_g_ = SD_USER_DEFINED;
    } else {
      casadi_error("Unknown linear solver for backward integration: " << iterative_solverB);
    }

    // Create a Jacobian if requested
    if (exact_jacobian_) {
      jac_ = getJac();
      alloc(jac_);
      alloc_w(jac_.sz_w() + jac_.nnz_out(0));
    }

    if (!jac_.isNull()) {
      casadi_assert_message(
                            jac_.size2_out(0)==jac_.size1_out(0),
                            "SundialsInterface::init: the jacobian of the forward problem must "
                            "be square but got " << jac_.output().dim());

      casadi_assert_message(
                            !jac_.sparsity_out(0).is_singular(),
                            "SundialsInterface::init: singularity - the jacobian of the forward "
                            "problem is structurally rank-deficient. sprank(J)="
                            << sprank(jac_.output()) << " (in stead of "<< jac_.size2_out(0)
                            << ")");
    }

    // Create a backwards Jacobian if requested
    if (exact_jacobianB_ && !g_.isNull()) jacB_ = getJacB();

    if (!jacB_.isNull()) {
      alloc(jacB_);
      alloc_w(jacB_.sz_w() + jacB_.nnz_out(0));
      casadi_assert_message(
                            jacB_.size2_out(0)==jacB_.size1_out(0),
                            "SundialsInterface::init: the jacobian of the backward problem must be "
                            "square but got " << jacB_.output().dim());

      casadi_assert_message(
                            !jacB_.sparsity_out(0).is_singular(),
                            "SundialsInterface::init: singularity - the jacobian of the backward"
                            " problem is structurally rank-deficient. sprank(J)="
                            << sprank(jacB_.output()) << " (instead of "
                            << jacB_.size2_out(0) << ")");
    }

    if (hasSetOption("linear_solver") && !jac_.isNull()) {
      // Options
      Dict linear_solver_options;
      if (hasSetOption("linear_solver_options")) {
        linear_solver_options = option("linear_solver_options");
      }

      // Create a linear solver
      linsol_ = Function::linsol("linsol", option("linear_solver"), jac_.sparsity_out(0),
                                 1, linear_solver_options);
      alloc(linsol_);
    }

    if ((hasSetOption("linear_solverB") || hasSetOption("linear_solver")) && !jacB_.isNull()) {
      // Linear solver options
      Dict opts;
      if (hasSetOption("linear_solver_optionsB")) {
        opts = option("linear_solver_optionsB");
      } else if (hasSetOption("linear_solver_options")) {
        opts = option("linear_solver_options");
      }

      // Create a linear solver
      std::string linear_solver_name =
        hasSetOption("linear_solverB") ? option("linear_solverB") : option("linear_solver");
      linsolB_ = Function::linsol("linsolB", linear_solver_name, jacB_.sparsity_out(0),
                                  1, opts);
      alloc(linsolB_);
    }

    // Allocate temporary memory
    //alloc_w(np_, true); // p_
    //alloc_w(nrp_, true); // rp_
  }

  void SundialsInterface::setup(Memory& m, const double**& arg, double**& res,
                                int*& iw, double*& w) {
    // Initialize the base classes
    Ivpsol::setup(m, arg, res, iw, w);

    // Get work vectors
    //p_ = w; w += np_;
    //rp_ = w; w += nrp_;
  }

  void SundialsInterface::reset(Memory& m, double t, const double* x,
                                const double* z, const double* p) {
    // Reset the base classes
    Ivpsol::reset(m, t, x, z, p);

    // Set the state
    casadi_copy(x, nx_, NV_DATA_S(xz_));
    casadi_copy(z, nz_, NV_DATA_S(xz_)+nx_);

    // Reset summation states
    N_VConst(0., q_);

    // Store parameters
    //casadi_copy(p, np_, p_);
  }

  void SundialsInterface::resetB(Memory& m, double t, const double* rx,
                                 const double* rz, const double* rp) {
    // Reset the base classes
    Ivpsol::resetB(m, t, rx, rz, rp);

    // Get the backward state
    casadi_copy(rx, nrx_, NV_DATA_S(rxz_));

    // Store parameters
    //casadi_copy(rp, nrp_, rp_);
  }

  std::pair<int, int> SundialsInterface::getBandwidth() const {
    std::pair<int, int> bw;

    // Get upper bandwidth
    if (hasSetOption("upper_bandwidth")) {
      bw.first = option("upper_bandwidth");
    } else if (!jac_.isNull()) {
      bw.first = jac_.getOutput().sparsity().bw_upper();
    } else {
      casadi_error("\"upper_bandwidth\" has not been set and cannot be "
                   "detected since exact Jacobian is not available.");
    }

    // Get lower bandwidth
    if (hasSetOption("lower_bandwidth")) {
      bw.second = option("lower_bandwidth");
    } else if (!jac_.isNull()) {
      bw.second = jac_.getOutput().sparsity().bw_lower();
    } else {
      casadi_error("\"lower_bandwidth\" has not been set and cannot be "
                   "detected since exact Jacobian is not available.");
    }

    return bw;
  }

  std::pair<int, int> SundialsInterface::getBandwidthB() const {
    std::pair<int, int> bw;

    // Get upper bandwidth
    if (hasSetOption("upper_bandwidthB")) {
      bw.first = option("upper_bandwidthB");
    } else if (!jacB_.isNull()) {
      bw.first = jacB_.getOutput().sparsity().bw_upper();
    } else {
      casadi_error("\"upper_bandwidthB\" has not been set and cannot be detected "
                   "since exact Jacobian for backward problem is not available.");
    }

    // Get lower bandwidth
    if (hasSetOption("lower_bandwidthB")) {
      bw.second = option("lower_bandwidthB");
    } else if (!jacB_.isNull()) {
      bw.second = jacB_.getOutput().sparsity().bw_lower();
    } else {
      casadi_error("\"lower_bandwidthB\" has not been set and cannot be detected "
                   "since exact Jacobian for backward problem is not available.");
    }

    return bw;
  }

} // namespace casadi


