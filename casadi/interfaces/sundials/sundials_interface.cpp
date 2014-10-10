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
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/mx/mx_tools.hpp"
#include "casadi/core/sx/sx_tools.hpp"
#include "casadi/core/function/mx_function.hpp"
#include "casadi/core/function/sx_function.hpp"

INPUTSCHEME(IntegratorInput)
OUTPUTSCHEME(IntegratorOutput)

using namespace std;
namespace casadi {

SundialsInterface::SundialsInterface(const Function& f, const Function& g)
    : IntegratorInternal(f, g) {
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
  addOption("linear_solver_options",       OT_DICTIONARY,       GenericType(),
            "Options to be passed to the linear solver");
  addOption("linear_solverB",              OT_STRING,     GenericType(),
            "A custom linear solver creator function for backwards integration "
            "[default: equal to linear_solver]");
  addOption("linear_solver_optionsB",      OT_DICTIONARY,       GenericType(),
            "Options to be passed to the linear solver for backwards integration "
            "[default: equal to linear_solver_options]");
}

SundialsInterface::~SundialsInterface() {
}

void SundialsInterface::init() {
  // Call the base class method
  IntegratorInternal::init();

  // Reset checkpoints counter
  ncheck_ = 0;

  // Read options
  abstol_ = getOption("abstol");
  reltol_ = getOption("reltol");
  exact_jacobian_ = getOption("exact_jacobian");
  exact_jacobianB_ = hasSetOption("exact_jacobianB") ?
      getOption("exact_jacobianB") && !g_.isNull() : exact_jacobian_;
  max_num_steps_ = getOption("max_num_steps");
  finite_difference_fsens_ = getOption("finite_difference_fsens");
  fsens_abstol_ =
      hasSetOption("fsens_abstol") ? static_cast<double>(getOption("fsens_abstol")) : abstol_;
  fsens_reltol_ =
      hasSetOption("fsens_reltol") ? static_cast<double>(getOption("fsens_reltol")) : reltol_;
  abstolB_ = hasSetOption("abstolB") ? static_cast<double>(getOption("abstolB")) : abstol_;
  reltolB_ = hasSetOption("reltolB") ? static_cast<double>(getOption("reltolB")) : reltol_;
  stop_at_end_ = getOption("stop_at_end");
  use_preconditioner_ = getOption("use_preconditioner");
  use_preconditionerB_ =  hasSetOption("use_preconditionerB") ?
      static_cast<bool>(getOption("use_preconditionerB")): use_preconditioner_;
  max_krylov_ = getOption("max_krylov");
  max_krylovB_ =
      hasSetOption("max_krylovB") ? static_cast<int>(getOption("max_krylovB")): max_krylov_;

  // Linear solver for forward integration
  if (getOption("linear_solver_type")=="dense") {
    linsol_f_ = SD_DENSE;
  } else if (getOption("linear_solver_type")=="banded") {
    linsol_f_ = SD_BANDED;
  } else if (getOption("linear_solver_type")=="iterative") {
    linsol_f_ = SD_ITERATIVE;

    // Iterative solver
    if (getOption("iterative_solver")=="gmres") {
      itsol_f_ = SD_GMRES;
    } else if (getOption("iterative_solver")=="bcgstab") {
      itsol_f_ = SD_BCGSTAB;
    } else if (getOption("iterative_solver")=="tfqmr") {
      itsol_f_ = SD_TFQMR;
    } else {
      throw CasadiException("Unknown sparse solver for forward integration");
    }

    // Preconditioning type
    if (getOption("pretype")=="none")               pretype_f_ = PREC_NONE;
    else if (getOption("pretype")=="left")          pretype_f_ = PREC_LEFT;
    else if (getOption("pretype")=="right")         pretype_f_ = PREC_RIGHT;
    else if (getOption("pretype")=="both")          pretype_f_ = PREC_BOTH;
    else
      throw CasadiException("Unknown preconditioning type for forward integration");
  } else if (getOption("linear_solver_type")=="user_defined") {
    linsol_f_ = SD_USER_DEFINED;
  } else {
    throw CasadiException("Unknown linear solver for forward integration");
  }


  std::string linear_solver_typeB = hasSetOption("linear_solver_typeB") ?
      getOption("linear_solver_typeB") : getOption("linear_solver_type");
  std::string iterative_solverB = hasSetOption("iterative_solverB") ?
      getOption("iterative_solverB") : getOption("iterative_solver");
  std::string pretypeB = hasSetOption("pretypeB") ? getOption("pretypeB"): getOption("pretype");

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
  if (exact_jacobian_) jac_ = getJac();
  // Initialize Jacobian if availabe
  if (!jac_.isNull() && !jac_.isInit()) jac_.init();

  if (!jac_.isNull()) {
    casadi_assert_message(
      jac_.output().size2()==jac_.output().size1(),
      "SundialsInterface::init: the jacobian of the forward problem must "
      "be square but got " << jac_.output().dimString());

    casadi_assert_message(
      !jac_.output().sparsity().isSingular(),
      "SundialsInterface::init: singularity - the jacobian of the forward "
      "problem is structurally rank-deficient. sprank(J)="
      << sprank(jac_.output()) << " (in stead of "<< jac_.output().size2()
      << ")");
  }

  // Create a backwards Jacobian if requested
  if (exact_jacobianB_ && !g_.isNull()) jacB_ = getJacB();
  // Initialize backwards  Jacobian if availabe
  if (!jacB_.isNull() && !jacB_.isInit()) jacB_.init();

  if (!jacB_.isNull()) {
    casadi_assert_message(
      jacB_.output().size2()==jacB_.output().size1(),
      "SundialsInterface::init: the jacobian of the backward problem must be "
      "square but got " << jacB_.output().dimString());

    casadi_assert_message(
      !jacB_.output().sparsity().isSingular(),
      "SundialsInterface::init: singularity - the jacobian of the backward"
      " problem is structurally rank-deficient. sprank(J)="
      << sprank(jacB_.output()) << " (instead of "
      << jacB_.output().size2() << ")");
  }

  if (hasSetOption("linear_solver") && !jac_.isNull()) {
    // Create a linear solver
    std::string linear_solver_name = getOption("linear_solver");
    linsol_ = LinearSolver(linear_solver_name, jac_.output().sparsity(), 1);
    // Pass options
    if (hasSetOption("linear_solver_options")) {
      linsol_.setOption(getOption("linear_solver_options"));
    }
    linsol_.init();
  }

  if ((hasSetOption("linear_solverB") || hasSetOption("linear_solver")) && !jacB_.isNull()) {
    // Create a linear solver
    std::string linear_solver_name =
        hasSetOption("linear_solverB") ? getOption("linear_solverB") : getOption("linear_solver");
    linsolB_ = LinearSolver(linear_solver_name, jacB_.output().sparsity(), 1);
    // Pass options
    if (hasSetOption("linear_solver_optionsB")) {
      linsolB_.setOption(getOption("linear_solver_optionsB"));
    } else if (hasSetOption("linear_solver_options")) {
      linsolB_.setOption(getOption("linear_solver_options"));
    }
    linsolB_.init();
  }
}

void SundialsInterface::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
  IntegratorInternal::deepCopyMembers(already_copied);
  linsol_ = deepcopy(linsol_, already_copied);
  linsolB_ = deepcopy(linsolB_, already_copied);
  jac_ = deepcopy(jac_, already_copied);
  jacB_ = deepcopy(jacB_, already_copied);
  f_fwd_ = deepcopy(f_fwd_, already_copied);
  g_fwd_ = deepcopy(g_fwd_, already_copied);
}

void SundialsInterface::reset() {
  // Reset the base classes
  IntegratorInternal::reset();
}

  std::pair<int, int> SundialsInterface::getBandwidth() const {
    std::pair<int, int> bw;

    // Get upper bandwidth
    if (hasSetOption("upper_bandwidth")) {
      bw.first = getOption("upper_bandwidth");
    } else if (!jac_.isNull()) {
      bw.first = jac_.getOutput().sparsity().bandwidthU();
    } else {
      casadi_error("\"upper_bandwidth\" has not been set and cannot be "
                   "detected since exact Jacobian is not available.");
    }

    // Get lower bandwidth
    if (hasSetOption("lower_bandwidth")) {
      bw.second = getOption("lower_bandwidth");
    } else if (!jac_.isNull()) {
      bw.second = jac_.getOutput().sparsity().bandwidthL();
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
      bw.first = getOption("upper_bandwidthB");
    } else if (!jacB_.isNull()) {
      bw.first = jacB_.getOutput().sparsity().bandwidthU();
    } else {
      casadi_error("\"upper_bandwidthB\" has not been set and cannot be detected "
                   "since exact Jacobian for backward problem is not available.");
    }

    // Get lower bandwidth
  if (hasSetOption("lower_bandwidthB")) {
      bw.second = getOption("lower_bandwidthB");
    } else if (!jacB_.isNull()) {
      bw.second = jacB_.getOutput().sparsity().bandwidthL();
    } else {
      casadi_error("\"lower_bandwidthB\" has not been set and cannot be detected "
                   "since exact Jacobian for backward problem is not available.");
    }

    return bw;
  }

} // namespace casadi


