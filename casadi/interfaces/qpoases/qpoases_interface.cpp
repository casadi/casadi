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


#include "qpoases_interface.hpp"

// Bug in qpOASES?
#define ALLOW_QPROBLEMB true
#define ALLOW_ALL_OPTIONS

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_QPSOL_QPOASES_EXPORT
  casadi_register_qpsol_qpoases(Qpsol::Plugin* plugin) {
    plugin->creator = QpoasesInterface::creator;
    plugin->name = "qpoases";
    plugin->doc = QpoasesInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_QPSOL_QPOASES_EXPORT casadi_load_qpsol_qpoases() {
    Qpsol::registerPlugin(casadi_register_qpsol_qpoases);
  }

  QpoasesInterface::QpoasesInterface(const std::string& name,
                                     const std::map<std::string, Sparsity>& st)
    : Qpsol(name, st) {

    addOption("nWSR", OT_INT,
              "The maximum number of working set recalculations to be performed during "
              "the initial homotopy. Default is 5(nx + nc)");
    addOption("CPUtime", OT_DOUBLE,
              "The maximum allowed CPU time in seconds for the whole initialisation"
              " (and the actually required one on output). Disabled if unset.");
    addOption("printLevel", OT_STRING,
              "Defines the amount of text output during QP solution, see Section 5.7");
    addOption("enableRamping", OT_BOOL,
              "Enables ramping.");
    addOption("enableFarBounds", OT_BOOL,
              "Enables the use of  far bounds.");
    addOption("enableFlippingBounds", OT_BOOL,
              "Enables the use of  flipping bounds.");
    addOption("enableRegularisation", OT_BOOL,
              "Enables automatic  Hessian regularisation.");
    addOption("enableFullLITests", OT_BOOL,
              "Enables condition-hardened  (but more expensive) LI test.");
    addOption("enableNZCTests", OT_BOOL,
              "Enables nonzero curvature  tests.");
    addOption("enableDriftCorrection", OT_INT,
              "Specifies the frequency of drift corrections: 0: turns them off.");
    addOption("enableCholeskyRefactorisation", OT_INT,
              "Specifies the frequency of a full re-factorisation of projected "
              "Hessian matrix: 0: turns them off,  1: uses them at each iteration etc.");
    addOption("enableEqualities", OT_BOOL,
              "Specifies whether equalities should be treated  as always active "
              "(True) or not (False)");
    addOption("terminationTolerance", OT_DOUBLE,
              "Relative termination tolerance to stop homotopy.");
    addOption("boundTolerance", OT_DOUBLE,
              "If upper and lower bounds differ less than this tolerance, they are regarded "
              "equal, i.e. as  equality constraint.");
    addOption("boundRelaxation", OT_DOUBLE,
              "Initial relaxation of bounds to start homotopy  and initial value for far bounds.");
    addOption("epsNum", OT_DOUBLE,
              "Numerator tolerance for ratio tests.");
    addOption("epsDen", OT_DOUBLE,
              "Denominator tolerance for ratio tests.");
    addOption("maxPrimalJump", OT_DOUBLE,
              "Maximum allowed jump in primal variables in  nonzero curvature tests.");
    addOption("maxDualJump", OT_DOUBLE,
              "Maximum allowed jump in dual variables in  linear independence tests.");
    addOption("initialRamping", OT_DOUBLE,
              "Start value for ramping strategy.");
    addOption("finalRamping", OT_DOUBLE,
              "Final value for ramping strategy.");
    addOption("initialFarBounds", OT_DOUBLE,
              "Initial size for far bounds.");
    addOption("growFarBounds", OT_DOUBLE,
              "Factor to grow far bounds.");
    addOption("initialStatusBounds", OT_STRING,
              "Initial status of bounds at first iteration.");
    addOption("epsFlipping", OT_DOUBLE,
              "Tolerance of squared Cholesky diagonal factor  which triggers flipping bound.");
    addOption("numRegularisationSteps", OT_INT,
              "Maximum number of successive regularisation steps.");
    addOption("epsRegularisation", OT_DOUBLE,
              "Scaling factor of identity matrix used for  Hessian regularisation.");
    addOption("numRefinementSteps", OT_INT,
              "Maximum number of iterative refinement steps.");
    addOption("epsIterRef", OT_DOUBLE,
              "Early termination tolerance for iterative  refinement.");
    addOption("epsLITests", OT_DOUBLE,
              "Tolerance for linear independence tests.");
    addOption("epsNZCTests", OT_DOUBLE,
              "Tolerance for nonzero curvature tests.");
  }

  QpoasesInterface::~QpoasesInterface() {
  }

  void QpoasesInterface::init(const Dict& opts) {
    Qpsol::init(opts);

    // Default options
    max_nWSR_ = 5 *(n_ + nc_);
    max_cputime_ = -1;
    ops_.setToDefault();

    // Read options
    for (auto&& op : opts) {
      if (op.first=="nWSR") {
        max_nWSR_ = op.second;
      } else if (op.first=="CPUtime") {
        max_cputime_ = op.second;
      } else if (op.first=="printLevel") {
        ops_.printLevel = to_PrintLevel(op.second);
      } else if (op.first=="enableRamping") {
        ops_.enableRamping = to_BooleanType(op.second);
      } else if (op.first=="enableFarBounds") {
        ops_.enableFarBounds = to_BooleanType(op.second);
      } else if (op.first=="enableFlippingBounds") {
        ops_.enableFlippingBounds = to_BooleanType(op.second);
      } else if (op.first=="enableRegularisation") {
        ops_.enableRegularisation = to_BooleanType(op.second);
      } else if (op.first=="enableFullLITests") {
        ops_.enableFullLITests = to_BooleanType(op.second);
      } else if (op.first=="enableNZCTests") {
        ops_.enableNZCTests = to_BooleanType(op.second);
      } else if (op.first=="enableDriftCorrection") {
        ops_.enableRegularisation = to_BooleanType(op.second);
      } else if (op.first=="enableCholeskyRefactorisation") {
        ops_.enableCholeskyRefactorisation = op.second;
      } else if (op.first=="terminationTolerance") {
        ops_.terminationTolerance = op.second;
      } else if (op.first=="boundTolerance") {
        ops_.boundTolerance = op.second;
      } else if (op.first=="boundRelaxation") {
        ops_.boundRelaxation = op.second;
      } else if (op.first=="epsNum") {
        ops_.epsNum = op.second;
      } else if (op.first=="epsDen") {
        ops_.epsDen = op.second;
      } else if (op.first=="maxPrimalJump") {
        ops_.maxPrimalJump = op.second;
      } else if (op.first=="maxDualJump") {
        ops_.maxDualJump = op.second;
      } else if (op.first=="initialRamping") {
        ops_.initialRamping = op.second;
      } else if (op.first=="finalRamping") {
        ops_.finalRamping = op.second;
      } else if (op.first=="initialFarBounds") {
        ops_.initialFarBounds = op.second;
      } else if (op.first=="growFarBounds") {
        ops_.growFarBounds = op.second;
      } else if (op.first=="initialStatusBounds") {
        ops_.initialStatusBounds = to_SubjectToStatus(op.second);
      } else if (op.first=="epsFlipping") {
        ops_.epsFlipping = op.second;
      } else if (op.first=="numRegularisationSteps") {
        ops_.numRegularisationSteps = op.second;
      } else if (op.first=="epsRegularisation") {
        ops_.epsRegularisation = op.second;
      } else if (op.first=="numRefinementSteps") {
        ops_.numRefinementSteps = op.second;
      } else if (op.first=="epsIterRef") {
        ops_.epsIterRef = op.second;
      } else if (op.first=="epsLITests") {
        ops_.epsLITests = op.second;
      } else if (op.first=="epsNZCTests") {
        ops_.epsNZCTests = op.second;
      }
    }

    // Allocate work vectors
    alloc_w(n_*n_, true); // h
    alloc_w(n_*nc_, true); // a
    alloc_w(n_, true); // g
    alloc_w(n_, true); // lbx
    alloc_w(n_, true); // ubx
    alloc_w(nc_, true); // lba
    alloc_w(nc_, true); // uba
    alloc_w(n_+nc_, true); // dual
  }

  void QpoasesInterface::init_memory(Memory& mem) const {
    QpoasesMemory& m = dynamic_cast<QpoasesMemory&>(mem);
    m.called_once = false;

    // Create qpOASES instance
    if (m.qp) delete m.qp;
    if (nc_==0) {
      m.qp = new qpOASES::QProblemB(n_);
    } else {
      m.qp = new qpOASES::SQProblem(n_, nc_);
    }

    // Pass to qpOASES
    m.qp->setOptions(ops_);
  }

  void QpoasesInterface::
  eval(Memory& mem, const double** arg, double** res, int* iw, double* w) const {
    QpoasesMemory& m = dynamic_cast<QpoasesMemory&>(mem);

    if (inputs_check_) {
      checkInputs(arg[QPSOL_LBX], arg[QPSOL_UBX], arg[QPSOL_LBA], arg[QPSOL_UBA]);
    }

    // Get quadratic term
    double* h=w; w += n_*n_;
    casadi_densify(arg[QPSOL_H], sparsity_in(QPSOL_H), h, false);

    // Get linear term
    double* a = w; w += n_*nc_;
    casadi_densify(arg[QPSOL_A], sparsity_in(QPSOL_A), a, true);

    // Maxiumum number of working set changes
    int nWSR = max_nWSR_;
    double cputime = max_cputime_;
    double *cputime_ptr = cputime<=0 ? 0 : &cputime;

    // Get the arguments to call qpOASES with
    double* g=w; w += n_;
    casadi_copy(arg[QPSOL_G], n_, g);
    double* lb=w; w += n_;
    casadi_copy(arg[QPSOL_LBX], n_, lb);
    double* ub=w; w += n_;
    casadi_copy(arg[QPSOL_UBX], n_, ub);
    double* lbA=w; w += nc_;
    casadi_copy(arg[QPSOL_LBA], nc_, lbA);
    double* ubA=w; w += nc_;
    casadi_copy(arg[QPSOL_UBA], nc_, ubA);

    int flag;
    if (!m.called_once) {
      if (nc_==0) {
        flag = static_cast<qpOASES::QProblemB*>(m.qp)->init(h, g, lb, ub, nWSR, cputime_ptr);
      } else {
        flag = static_cast<qpOASES::SQProblem*>(m.qp)->init(h, g, a, lb, ub, lbA, ubA,
                                                           nWSR, cputime_ptr);
      }
      m.called_once = true;
    } else {
      if (nc_==0) {
        static_cast<qpOASES::QProblemB*>(m.qp)->reset();
        flag = static_cast<qpOASES::QProblemB*>(m.qp)->init(h, g, lb, ub, nWSR, cputime_ptr);
        //flag = static_cast<qpOASES::QProblemB*>(m.qp)->hotstart(g, lb, ub, nWSR, cputime_ptr);
      } else {
        flag = static_cast<qpOASES::SQProblem*>(m.qp)->hotstart(h, g, a, lb, ub, lbA, ubA,
                                                               nWSR, cputime_ptr);
      }
    }
    if (flag!=qpOASES::SUCCESSFUL_RETURN && flag!=qpOASES::RET_MAX_NWSR_REACHED) {
      throw CasadiException("qpOASES failed: " + getErrorMessage(flag));
    }

    // Get optimal cost
    if (res[QPSOL_COST]) *res[QPSOL_COST] = m.qp->getObjVal();

    // Get the primal solution
    if (res[QPSOL_X]) m.qp->getPrimalSolution(res[QPSOL_X]);

    // Get the dual solution
    if (res[QPSOL_LAM_X] || res[QPSOL_LAM_A]) {
      double* dual=w; w += n_+nc_;
      m.qp->getDualSolution(dual);
      casadi_scal(n_+nc_, -1., dual);
      casadi_copy(dual, n_, res[QPSOL_LAM_X]);
      casadi_copy(dual+n_, nc_, res[QPSOL_LAM_A]);
    }
  }

  std::string QpoasesInterface::getErrorMessage(int flag) {
    switch (flag) {
    case qpOASES::SUCCESSFUL_RETURN:
      return "Successful return.";
    case qpOASES::RET_DIV_BY_ZERO:
      return "Division by zero.";
    case qpOASES::RET_INDEX_OUT_OF_BOUNDS:
      return "Index out of bounds.";
    case qpOASES::RET_INVALID_ARGUMENTS:
      return "At least one of the arguments is invalid.";
    case qpOASES::RET_ERROR_UNDEFINED:
      return "Error number undefined.";
    case qpOASES::RET_WARNING_UNDEFINED:
      return "Warning number undefined.";
    case qpOASES::RET_INFO_UNDEFINED:
      return "Info number undefined.";
    case qpOASES::RET_EWI_UNDEFINED:
      return "Error/warning/info number undefined.";
    case qpOASES::RET_AVAILABLE_WITH_LINUX_ONLY:
      return "This function is available under Linux only.";
    case qpOASES::RET_UNKNOWN_BUG:
      return "The error occured is not yet known.";
    case qpOASES::RET_PRINTLEVEL_CHANGED:
      return "Print level changed.";
    case qpOASES::RET_NOT_YET_IMPLEMENTED:
      return "Requested function is not yet implemented in this version of qpOASES.";
      // Indexlist
    case qpOASES::RET_INDEXLIST_MUST_BE_REORDERD:
      return "Index list has to be reordered.";
    case qpOASES::RET_INDEXLIST_EXCEEDS_MAX_LENGTH:
      return "Index list exceeds its maximal physical length.";
    case qpOASES::RET_INDEXLIST_CORRUPTED:
      return "Index list corrupted.";
    case qpOASES::RET_INDEXLIST_OUTOFBOUNDS:
      return "Physical index is out of bounds.";
    case qpOASES::RET_INDEXLIST_ADD_FAILED:
      return "Adding indices from another index set failed.";
    case qpOASES::RET_INDEXLIST_INTERSECT_FAILED:
      return "Intersection with another index set failed.";
      // SubjectTo / Bounds / Constraints
    case qpOASES::RET_INDEX_ALREADY_OF_DESIRED_STATUS:
      return "Index is already of desired status.";
    case qpOASES::RET_ADDINDEX_FAILED:
      return "Adding index to index set failed.";
    case qpOASES::RET_REMOVEINDEX_FAILED:
      return "Removing index from index set failed.";
    case qpOASES::RET_SWAPINDEX_FAILED:
      return "Cannot swap between different indexsets.";
    case qpOASES::RET_NOTHING_TO_DO:
      return "Nothing to do.";
    case qpOASES::RET_SETUP_BOUND_FAILED:
      return "Setting up bound index failed.";
    case qpOASES::RET_SETUP_CONSTRAINT_FAILED:
      return "Setting up constraint index failed.";
    case qpOASES::RET_MOVING_BOUND_FAILED:
      return "Moving bound between index sets failed.";
    case qpOASES::RET_MOVING_CONSTRAINT_FAILED:
      return "Moving constraint between index sets failed.";
    case qpOASES::RET_SHIFTING_FAILED:
      return "Shifting of bounds/constraints failed.";
    case qpOASES::RET_ROTATING_FAILED:
      return "Rotating of bounds/constraints failed.";
      // QProblem
    case qpOASES::RET_QPOBJECT_NOT_SETUP:
      return "The QP object has not been setup correctly, use another constructor.";
    case qpOASES::RET_QP_ALREADY_INITIALISED:
      return "QProblem has already been initialized.";
    case qpOASES::RET_NO_INIT_WITH_STANDARD_SOLVER:
      return "Initialisation via extern QP solver is not yet implemented.";
    case qpOASES::RET_RESET_FAILED:
      return "Reset failed.";
    case qpOASES::RET_INIT_FAILED:
      return "Initialisation failed.";
    case qpOASES::RET_INIT_FAILED_TQ:
      return "Initialisation failed due to TQ factorisation.";
    case qpOASES::RET_INIT_FAILED_CHOLESKY:
      return "Initialisation failed due to Cholesky decomposition.";
    case qpOASES::RET_INIT_FAILED_HOTSTART:
      return "Initialisation failed! QP could not be solved!";
    case qpOASES::RET_INIT_FAILED_INFEASIBILITY:
      return "Initial QP could not be solved due to infeasibility!";
    case qpOASES::RET_INIT_FAILED_UNBOUNDEDNESS:
      return "Initial QP could not be solved due to unboundedness!";
    case qpOASES::RET_INIT_SUCCESSFUL:
      return "Initialisation done.";
    case qpOASES::RET_OBTAINING_WORKINGSET_FAILED:
      return "Failed to obtain working set for auxiliary QP.";
    case qpOASES::RET_SETUP_WORKINGSET_FAILED:
      return "Failed to setup working set for auxiliary QP.";
    case qpOASES::RET_SETUP_AUXILIARYQP_FAILED:
      return "Failed to setup auxiliary QP for initialized homotopy.";
    case qpOASES::RET_NO_EXTERN_SOLVER:
      return "No extern QP solver available.";
    case qpOASES::RET_QP_UNBOUNDED:
      return "QP is unbounded.";
    case qpOASES::RET_QP_INFEASIBLE:
      return "QP is infeasible.";
    case qpOASES::RET_QP_NOT_SOLVED:
      return "Problems occured while solving QP with standard solver.";
    case qpOASES::RET_QP_SOLVED:
      return "QP successfully solved.";
    case qpOASES::RET_UNABLE_TO_SOLVE_QP:
      return "Problems occured while solving QP.";
    case qpOASES::RET_INITIALISATION_STARTED:
      return "Starting problem initialisation.";
    case qpOASES::RET_HOTSTART_FAILED:
      return "Unable to perform homotopy due to internal error.";
    case qpOASES::RET_HOTSTART_FAILED_TO_INIT:
      return "Unable to initialise problem.";
    case qpOASES::RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED:
      return "Unable to perform homotopy as previous QP is not solved.";
    case qpOASES::RET_ITERATION_STARTED:
      return "Iteration...";
    case qpOASES::RET_SHIFT_DETERMINATION_FAILED:
      return "Determination of shift of the QP data failed.";
    case qpOASES::RET_STEPDIRECTION_DETERMINATION_FAILED:
      return "Determination of step direction failed.";
    case qpOASES::RET_STEPLENGTH_DETERMINATION_FAILED:
      return "Determination of step direction failed.";
    case qpOASES::RET_OPTIMAL_SOLUTION_FOUND:
      return "Optimal solution of neighbouring QP found.";
    case qpOASES::RET_HOMOTOPY_STEP_FAILED:
      return "Unable to perform homotopy step.";
    case qpOASES::RET_HOTSTART_STOPPED_INFEASIBILITY:
      return "Premature homotopy termination because QP is infeasible.";
    case qpOASES::RET_HOTSTART_STOPPED_UNBOUNDEDNESS:
      return "Premature homotopy termination because QP is unbounded.";
    case qpOASES::RET_WORKINGSET_UPDATE_FAILED:
      return "Unable to update working sets according to initial guesses.";
    case qpOASES::RET_MAX_NWSR_REACHED:
      return "Maximum number of working set recalculations performed.";
    case qpOASES::RET_CONSTRAINTS_NOT_SPECIFIED:
      return "Problem does comprise constraints! "
        "You also have to specify new constraints' bounds.";
    case qpOASES::RET_INVALID_FACTORISATION_FLAG:
      return "Invalid factorisation flag.";
    case qpOASES::RET_UNABLE_TO_SAVE_QPDATA:
      return "Unable to save QP data.";
    case qpOASES::RET_STEPDIRECTION_FAILED_TQ:
      return "Abnormal termination due to TQ factorisation.";
    case qpOASES::RET_STEPDIRECTION_FAILED_CHOLESKY:
      return "Abnormal termination due to Cholesky factorisation.";
    case qpOASES::RET_CYCLING_DETECTED:
      return "Cycling detected.";
    case qpOASES::RET_CYCLING_NOT_RESOLVED:
      return "Cycling cannot be resolved, QP probably infeasible.";
    case qpOASES::RET_CYCLING_RESOLVED:
      return "Cycling probably resolved.";
    case qpOASES::RET_STEPSIZE:
      return "For displaying performed stepsize.";
    case qpOASES::RET_STEPSIZE_NONPOSITIVE:
      return "For displaying non-positive stepsize.";
    case qpOASES::RET_SETUPSUBJECTTOTYPE_FAILED:
      return "Setup of SubjectToTypes failed.";
    case qpOASES::RET_ADDCONSTRAINT_FAILED:
      return "Addition of constraint to working set failed.";
    case qpOASES::RET_ADDCONSTRAINT_FAILED_INFEASIBILITY:
      return "Addition of constraint to working set failed (due to QP infeasibility).";
    case qpOASES::RET_ADDBOUND_FAILED:
      return "Addition of bound to working set failed.";
    case qpOASES::RET_ADDBOUND_FAILED_INFEASIBILITY:
      return "Addition of bound to working set failed (due to QP infeasibility).";
    case qpOASES::RET_REMOVECONSTRAINT_FAILED:
      return "Removal of constraint from working set failed.";
    case qpOASES::RET_REMOVEBOUND_FAILED:
      return "Removal of bound from working set failed.";
    case qpOASES::RET_REMOVE_FROM_ACTIVESET:
      return "Removing from active set...";
    case qpOASES::RET_ADD_TO_ACTIVESET:
      return "Adding to active set...";
    case qpOASES::RET_REMOVE_FROM_ACTIVESET_FAILED:
      return "Removing from active set failed.";
    case qpOASES::RET_ADD_TO_ACTIVESET_FAILED:
      return "Adding to active set failed.";
    case qpOASES::RET_CONSTRAINT_ALREADY_ACTIVE:
      return "Constraint is already active.";
    case qpOASES::RET_ALL_CONSTRAINTS_ACTIVE:
      return "All constraints are active, no further constraint can be added.";
    case qpOASES::RET_LINEARLY_DEPENDENT:
      return "New bound/constraint is linearly dependent.";
    case qpOASES::RET_LINEARLY_INDEPENDENT:
      return "New bound/constraint is linearly independent.";
    case qpOASES::RET_LI_RESOLVED:
      return "Linear indepence of active contraint matrix successfully resolved.";
    case qpOASES::RET_ENSURELI_FAILED:
      return "Failed to ensure linear indepence of active contraint matrix.";
    case qpOASES::RET_ENSURELI_FAILED_TQ:
      return "Abnormal termination due to TQ factorisation.";
    case qpOASES::RET_ENSURELI_FAILED_NOINDEX:
      return "QP is infeasible.";
    case qpOASES::RET_ENSURELI_FAILED_CYCLING:
      return "QP is infeasible.";
    case qpOASES::RET_BOUND_ALREADY_ACTIVE:
      return "Bound is already active.";
    case qpOASES::RET_ALL_BOUNDS_ACTIVE:
      return "All bounds are active, no further bound can be added.";
    case qpOASES::RET_CONSTRAINT_NOT_ACTIVE:
      return "Constraint is not active.";
    case qpOASES::RET_BOUND_NOT_ACTIVE:
      return "Bound is not active.";
    case qpOASES::RET_HESSIAN_NOT_SPD:
      return "Projected Hessian matrix not positive definite.";
    case qpOASES::RET_HESSIAN_INDEFINITE:
      return "Hessian matrix is indefinite.";
    case qpOASES::RET_MATRIX_SHIFT_FAILED:
      return "Unable to update matrices or to transform vectors.";
    case qpOASES::RET_MATRIX_FACTORISATION_FAILED:
      return "Unable to calculate new matrix factorisations.";
    case qpOASES::RET_PRINT_ITERATION_FAILED:
      return "Unable to print information on current iteration.";
    case qpOASES::RET_NO_GLOBAL_MESSAGE_OUTPUTFILE:
      return "No global message output file initialized.";
    case qpOASES::RET_DISABLECONSTRAINTS_FAILED:
      return "Unable to disbable constraints.";
    case qpOASES::RET_ENABLECONSTRAINTS_FAILED:
      return "Unable to enbable constraints.";
    case qpOASES::RET_ALREADY_ENABLED:
      return "Bound or constraint is already enabled.";
    case qpOASES::RET_ALREADY_DISABLED:
      return "Bound or constraint is already disabled.";
    case qpOASES::RET_NO_HESSIAN_SPECIFIED:
      return "No Hessian matrix has been specified.";
    case qpOASES::RET_USING_REGULARISATION:
      return "Using regularisation as Hessian matrix is not positive definite.";
    case qpOASES::RET_EPS_MUST_BE_POSITVE:
      return "Eps for regularisation must be sufficiently positive.";
    case qpOASES::RET_REGSTEPS_MUST_BE_POSITVE:
      return "Maximum number of regularisation steps must be non-negative.";
    case qpOASES::RET_HESSIAN_ALREADY_REGULARISED:
      return "Hessian has been already regularised.";
    case qpOASES::RET_CANNOT_REGULARISE_IDENTITY:
      return "Identity Hessian matrix cannot be regularised.";
    case qpOASES::RET_NO_REGSTEP_NWSR:
      return "No additional regularisation step could be performed due to limits.";
    case qpOASES::RET_FEWER_REGSTEPS_NWSR:
      return "Fewer additional regularisation steps have been performed due to limits.";
    case qpOASES::RET_CHOLESKY_OF_ZERO_HESSIAN:
      return "Cholesky decomposition of (unregularised) zero Hessian matrix.";
    case qpOASES::RET_CONSTRAINTS_ARE_NOT_SCALED:
      return "When defining __MANY_CONSTRAINTS__, l1 norm of each "
        "constraint must be not greater than one.";
    case qpOASES::RET_ERROR_IN_CONSTRAINTPRODUCT:
      return "Error in user-defined constraint product function.";
      // SQProblem
    case qpOASES::RET_UPDATEMATRICES_FAILED:
      return "Unable to update QP matrices.";
    case qpOASES::RET_UPDATEMATRICES_FAILED_AS_QP_NOT_SOLVED:
      return "Unable to update matrices as previous QP is not solved.";
      // Utils
    case qpOASES::RET_UNABLE_TO_OPEN_FILE:
      return "Unable to open file.";
    case qpOASES::RET_UNABLE_TO_WRITE_FILE:
      return "Unable to write into file.";
    case qpOASES::RET_UNABLE_TO_READ_FILE:
      return "Unable to read from file.";
    case qpOASES::RET_FILEDATA_INCONSISTENT:
      return "File contains inconsistent data.";
      // SolutionAnalysis
    case qpOASES::RET_UNABLE_TO_ANALYSE_QPROBLEM:
      return "Unable to analyse (S)QProblem(B) object";
      // Benchmark
    case qpOASES::RET_NWSR_SET_TO_ONE:
      return "Maximum number of working set changes was set to 1.";
    case qpOASES::RET_BENCHMARK_ABORTED:
      return "Benchmark aborted.";
    case qpOASES::RET_UNABLE_TO_READ_BENCHMARK:
      return "Unable to read benchmark data.";
    case qpOASES::RET_INITIAL_QP_SOLVED:
      return "Initial QP solved.";
    case qpOASES::RET_QP_SOLUTION_STARTED:
      return "Solving QP...";
    case qpOASES::RET_BENCHMARK_SUCCESSFUL:
      return "Benchmark terminated successfully.";
    }

    // Default error message
    stringstream ss;
    ss << "Unknown error flag: " << flag << ". Consult qpOASES documentation.";
    return ss.str();
  }

  bool QpoasesInterface::from_BooleanType(qpOASES::BooleanType b) {
    switch (b) {
    case qpOASES::BT_TRUE:              return true;
    case qpOASES::BT_FALSE:             return false;
    }
    casadi_error("not_implemented");
  }

  qpOASES::BooleanType QpoasesInterface::to_BooleanType(bool b) {
    return b ? qpOASES::BT_TRUE : qpOASES::BT_FALSE;
  }

  std::string QpoasesInterface::from_SubjectToStatus(qpOASES::SubjectToStatus b) {
    switch (b) {
    case qpOASES::ST_INACTIVE:          return "inactive";
    case qpOASES::ST_LOWER:             return "lower";
    case qpOASES::ST_UPPER:             return "upper";
    case qpOASES::ST_INFEASIBLE_LOWER:  return "infeasible_lower";
    case qpOASES::ST_INFEASIBLE_UPPER:  return "infeasible_upper";
    case qpOASES::ST_UNDEFINED:         return "undefined";
    }
    casadi_error("not_implemented");
  }

  qpOASES::SubjectToStatus QpoasesInterface::to_SubjectToStatus(std::string b) {
    if (b == "inactive") {
      return qpOASES::ST_INACTIVE;
    } else if (b == "lower") {
      return qpOASES::ST_LOWER;
    } else if (b == "infeasible_lower") {
      return qpOASES::ST_INFEASIBLE_LOWER;
    } else if (b == "infeasible_upper") {
      return qpOASES::ST_INFEASIBLE_UPPER;
    } else if (b == "undefined") {
      return qpOASES::ST_UNDEFINED;
    } else {
      casadi_error("No such qpOASES::SubjectToStatus: " + b);
    }
  }

  std::string QpoasesInterface::from_PrintLevel(qpOASES::PrintLevel b) {
    switch (b) {
    case qpOASES::PL_TABULAR:           return "tabular";
    case qpOASES::PL_NONE:              return "none";
    case qpOASES::PL_LOW:               return "low";
    case qpOASES::PL_MEDIUM:            return "medium";
    case qpOASES::PL_HIGH:              return "high";
    }
    casadi_error("not_implemented");
  }

  qpOASES::PrintLevel QpoasesInterface::to_PrintLevel(std::string b) {
    if (b == "tabular") {
      return qpOASES::PL_TABULAR;
    } else if (b == "none") {
      return qpOASES::PL_NONE;
    } else if (b == "low") {
      return qpOASES::PL_LOW;
    } else if (b == "medium") {
      return qpOASES::PL_MEDIUM;
    } else if (b == "high") {
      return qpOASES::PL_HIGH;
    } else {
      casadi_error("No such qpOASES::PrintLevel: " + b);
    }
  }

  QpoasesMemory::QpoasesMemory() {
    this->qp = 0;
  }

  QpoasesMemory::~QpoasesMemory() {
    if (this->qp) delete this->qp;
  }

} // namespace casadi
