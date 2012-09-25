/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include "qpoases_internal.hpp"
#include "../../casadi/fx/qp_solver.hpp"

#include "../../casadi/stl_vector_tools.hpp"
#include "../../casadi/matrix/matrix_tools.hpp"

const bool ALLOW_QPROBLEMB = true;

using namespace std;
namespace CasADi {
namespace Interfaces {

QPOasesInternal* QPOasesInternal::clone() const{
  // Return a deep copy
  QPOasesInternal* node = new QPOasesInternal(input(QP_H).sparsity(),input(QP_A).sparsity());
  if(!node->is_init_)
    node->init();
  return node;
}
  
QPOasesInternal::QPOasesInternal(const CRSSparsity& H, const CRSSparsity& A) : QPSolverInternal(H,A){
  addOption("nWSR",       OT_INTEGER,     GenericType(), "The maximum number of working set recalculations to be performed during the initial homotopy. Default is 5(nx + nc)");
  addOption("CPUtime",    OT_REAL,        GenericType(), "The maximum allowed CPU time in seconds for the whole initialisation (and the actually required one on output). Disabled if unset.");
  addOption("printLevel", OT_STRING,      GenericType(), "Defines the amount of text output during QP solution, see Section 5.7","none|low|medium|high");
  
  called_once_ = false;
  qp_ = 0;
}

QPOasesInternal::~QPOasesInternal(){ 
  if(qp_!=0) delete qp_;
}

void QPOasesInternal::init(){
  QPSolverInternal::init();

  // Read options
  if(hasSetOption("nWSR")){
    max_nWSR_ = getOption("nWSR");
    casadi_assert(max_nWSR_>=0);
  } else {
    max_nWSR_ = 5 *(nx_ + nc_);
  }

  if(hasSetOption("CPUtime")){
    max_cputime_ = getOption("CPUtime");
    casadi_assert(max_cputime_>0);
  } else {
    max_cputime_ = -1;
  }
  
  // Create data for H if not dense
  if(!input(QP_H).sparsity().dense()) h_data_.resize(nx_*nx_);
  
  // Create data for A if not dense
  if(!input(QP_A).sparsity().dense()) a_data_.resize(nx_*nc_);
  
  // Dual solution vector
  dual_.resize(nx_+nc_);
  
  // Create qpOASES instance
  if(qp_) delete qp_;
  if(ALLOW_QPROBLEMB && nc_==0){
    qp_ = new qpOASES::QProblemB(nx_);
  } else {
    qp_ = new qpOASES::SQProblem(nx_,nc_);
  }
  called_once_ = false;

  // Set options
  //qpOASES::Options myoptions;
  //myoptions.setToDefault(); 
  
  // Get print level
  if(hasSetOption("printLevel")){
    if(getOption("printLevel")=="none"){
      qp_->setPrintLevel(qpOASES::PL_NONE);
      //myoptions.printLevel = qpOASES::PL_NONE;
    } else if(getOption("printLevel")=="low"){
      qp_->setPrintLevel(qpOASES::PL_LOW);
      //myoptions.printLevel = qpOASES::PL_LOW;
    } else if(getOption("printLevel")=="medium"){
      qp_->setPrintLevel(qpOASES::PL_MEDIUM);
      //myoptions.printLevel = qpOASES::PL_MEDIUM;
    } else if(getOption("printLevel")=="high"){
      qp_->setPrintLevel(qpOASES::PL_HIGH);
      //myoptions.printLevel = qpOASES::PL_HIGH;
    }
  }
  
  // Pass to qpOASES
  //qp_->setOption(myoptions);
}

void QPOasesInternal::evaluate(int nfdir, int nadir) {
  if (nfdir!=0 || nadir!=0) throw CasadiException("QPOasesInternal::evaluate() not implemented for forward or backward mode");

  if(verbose()){
//     cout << "X_INIT = " << input(QP_X_INIT) << endl;
//     cout << "LAMBDA_INIT = " << input(QP_LAMBDA_INIT) << endl;
    cout << "LBX = " << input(QP_LBX) << endl;
    cout << "UBX = " << input(QP_UBX) << endl;
    cout << "LBA = " << input(QP_LBA) << endl;
    cout << "UBA = " << input(QP_UBA) << endl;
  }
  
  // Get pointer to H
  const double* h=0;
  if(h_data_.empty()){
    // No copying needed
    h = getPtr(input(QP_H));
  } else {
    // First copy to dense array
    input(QP_H).get(h_data_,DENSE);
    h = getPtr(h_data_);
  }
  
  // Get pointer to A
  const double* a=0;
  if(a_data_.empty()){
    // No copying needed
    a = getPtr(input(QP_A));
  } else {
    // First copy to dense array
    input(QP_A).get(a_data_,DENSE);
    a = getPtr(a_data_);
  }
  
  // Maxiumum number of working set changes
  int nWSR = max_nWSR_;
  double cputime = max_cputime_;
  double *cputime_ptr = cputime<=0 ? 0 : &cputime;

  // Get the arguments to call qpOASES with
  const double* g = getPtr(input(QP_G));
  const double* lb = getPtr(input(QP_LBX));
  const double* ub = getPtr(input(QP_UBX));
  const double* lbA = getPtr(input(QP_LBA));
  const double* ubA = getPtr(input(QP_UBA));

  int flag;
  if(!called_once_){
    if(ALLOW_QPROBLEMB && nc_==0){
      flag = static_cast<qpOASES::QProblemB*>(qp_)->init(h,g,lb,ub,nWSR,cputime_ptr);
    } else {
      flag = static_cast<qpOASES::SQProblem*>(qp_)->init(h,g,a,lb,ub,lbA,ubA,nWSR,cputime_ptr);
    }
    called_once_ = true;
  } else {
    if(ALLOW_QPROBLEMB && nc_==0){
      static_cast<qpOASES::QProblemB*>(qp_)->reset();
      flag = static_cast<qpOASES::QProblemB*>(qp_)->init(h,g,lb,ub,nWSR,cputime_ptr);
      //flag = static_cast<qpOASES::QProblemB*>(qp_)->hotstart(g,lb,ub,nWSR,cputime_ptr);
    } else {
      flag = static_cast<qpOASES::SQProblem*>(qp_)->hotstart(h,g,a,lb,ub,lbA,ubA,nWSR, cputime_ptr);
    }
  }
  if(flag!=qpOASES::SUCCESSFUL_RETURN && flag!=qpOASES::RET_MAX_NWSR_REACHED){
    throw CasadiException("qpOASES failed: " + getErrorMessage(flag));
  }

  // Get optimal cost
  output(QP_COST).set(qp_->getObjVal());

  // Get the primal solution
  qp_->getPrimalSolution(&output(QP_PRIMAL).front());
  
  // Get the dual solution
  qp_->getDualSolution(&dual_.front());
  
  // Split up the dual solution in multipliers for the simple bounds and the linear bounds
  transform(dual_.begin(),   dual_.begin()+nx_,output(QP_LAMBDA_X).begin(),negate<double>());
  transform(dual_.begin()+nx_,dual_.end(),     output(QP_LAMBDA_A).begin(),negate<double>());
}

std::string QPOasesInternal::getErrorMessage(int flag){
  switch(flag){
    case qpOASES::SUCCESSFUL_RETURN: return "Successful return.";
    case qpOASES::RET_DIV_BY_ZERO: return "Division by zero.";
    case qpOASES::RET_INDEX_OUT_OF_BOUNDS: return "Index out of bounds.";
    case qpOASES::RET_INVALID_ARGUMENTS: return "At least one of the arguments is invalid.";
    case qpOASES::RET_ERROR_UNDEFINED: return "Error number undefined.";
    case qpOASES::RET_WARNING_UNDEFINED: return "Warning number undefined.";
    case qpOASES::RET_INFO_UNDEFINED: return "Info number undefined.";
    case qpOASES::RET_EWI_UNDEFINED: return "Error/warning/info number undefined.";
    case qpOASES::RET_AVAILABLE_WITH_LINUX_ONLY: return "This function is available under Linux only.";
    case qpOASES::RET_UNKNOWN_BUG: return "The error occured is not yet known.";
    case qpOASES::RET_PRINTLEVEL_CHANGED: return "Print level changed.";
    case qpOASES::RET_NOT_YET_IMPLEMENTED: return "Requested function is not yet implemented in this version of qpOASES.";
    // Indexlist
    case qpOASES::RET_INDEXLIST_MUST_BE_REORDERD: return "Index list has to be reordered.";
    case qpOASES::RET_INDEXLIST_EXCEEDS_MAX_LENGTH: return "Index list exceeds its maximal physical length.";
    case qpOASES::RET_INDEXLIST_CORRUPTED: return "Index list corrupted.";
    case qpOASES::RET_INDEXLIST_OUTOFBOUNDS: return "Physical index is out of bounds.";
    case qpOASES::RET_INDEXLIST_ADD_FAILED: return "Adding indices from another index set failed.";
    case qpOASES::RET_INDEXLIST_INTERSECT_FAILED: return "Intersection with another index set failed.";
    // SubjectTo / Bounds / Constraints
    case qpOASES::RET_INDEX_ALREADY_OF_DESIRED_STATUS: return "Index is already of desired status.";
    case qpOASES::RET_ADDINDEX_FAILED: return "Adding index to index set failed.";
    case qpOASES::RET_REMOVEINDEX_FAILED: return "Removing index from index set failed.";
    case qpOASES::RET_SWAPINDEX_FAILED: return "Cannot swap between different indexsets.";
    case qpOASES::RET_NOTHING_TO_DO: return "Nothing to do.";
    case qpOASES::RET_SETUP_BOUND_FAILED: return "Setting up bound index failed.";
    case qpOASES::RET_SETUP_CONSTRAINT_FAILED: return "Setting up constraint index failed.";
    case qpOASES::RET_MOVING_BOUND_FAILED: return "Moving bound between index sets failed.";
    case qpOASES::RET_MOVING_CONSTRAINT_FAILED: return "Moving constraint between index sets failed.";
    case qpOASES::RET_SHIFTING_FAILED: return "Shifting of bounds/constraints failed.";
    case qpOASES::RET_ROTATING_FAILED: return "Rotating of bounds/constraints failed.";
    // QProblem
    case qpOASES::RET_QPOBJECT_NOT_SETUP: return "The QP object has not been setup correctly, use another constructor.";
    case qpOASES::RET_QP_ALREADY_INITIALISED: return "QProblem has already been initialised.";
    case qpOASES::RET_NO_INIT_WITH_STANDARD_SOLVER: return "Initialisation via extern QP solver is not yet implemented.";
    case qpOASES::RET_RESET_FAILED: return "Reset failed.";
    case qpOASES::RET_INIT_FAILED: return "Initialisation failed.";
    case qpOASES::RET_INIT_FAILED_TQ: return "Initialisation failed due to TQ factorisation.";
    case qpOASES::RET_INIT_FAILED_CHOLESKY: return "Initialisation failed due to Cholesky decomposition.";
    case qpOASES::RET_INIT_FAILED_HOTSTART: return "Initialisation failed! QP could not be solved!";
    case qpOASES::RET_INIT_FAILED_INFEASIBILITY: return "Initial QP could not be solved due to infeasibility!";
    case qpOASES::RET_INIT_FAILED_UNBOUNDEDNESS: return "Initial QP could not be solved due to unboundedness!";
    case qpOASES::RET_INIT_SUCCESSFUL: return "Initialisation done.";
    case qpOASES::RET_OBTAINING_WORKINGSET_FAILED: return "Failed to obtain working set for auxiliary QP.";
    case qpOASES::RET_SETUP_WORKINGSET_FAILED: return "Failed to setup working set for auxiliary QP.";
    case qpOASES::RET_SETUP_AUXILIARYQP_FAILED: return "Failed to setup auxiliary QP for initialised homotopy.";
    case qpOASES::RET_NO_EXTERN_SOLVER: return "No extern QP solver available.";
    case qpOASES::RET_QP_UNBOUNDED: return "QP is unbounded.";
    case qpOASES::RET_QP_INFEASIBLE: return "QP is infeasible.";
    case qpOASES::RET_QP_NOT_SOLVED: return "Problems occured while solving QP with standard solver.";
    case qpOASES::RET_QP_SOLVED: return "QP successfully solved.";
    case qpOASES::RET_UNABLE_TO_SOLVE_QP: return "Problems occured while solving QP.";
    case qpOASES::RET_INITIALISATION_STARTED: return "Starting problem initialisation.";
    case qpOASES::RET_HOTSTART_FAILED: return "Unable to perform homotopy due to internal error.";
    case qpOASES::RET_HOTSTART_FAILED_TO_INIT: return "Unable to initialise problem.";
    case qpOASES::RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED: return "Unable to perform homotopy as previous QP is not solved.";
    case qpOASES::RET_ITERATION_STARTED: return "Iteration...";
    case qpOASES::RET_SHIFT_DETERMINATION_FAILED: return "Determination of shift of the QP data failed.";
    case qpOASES::RET_STEPDIRECTION_DETERMINATION_FAILED: return "Determination of step direction failed.";
    case qpOASES::RET_STEPLENGTH_DETERMINATION_FAILED: return "Determination of step direction failed.";
    case qpOASES::RET_OPTIMAL_SOLUTION_FOUND: return "Optimal solution of neighbouring QP found.";
    case qpOASES::RET_HOMOTOPY_STEP_FAILED: return "Unable to perform homotopy step.";
    case qpOASES::RET_HOTSTART_STOPPED_INFEASIBILITY: return "Premature homotopy termination because QP is infeasible.";
    case qpOASES::RET_HOTSTART_STOPPED_UNBOUNDEDNESS: return "Premature homotopy termination because QP is unbounded.";
    case qpOASES::RET_WORKINGSET_UPDATE_FAILED: return "Unable to update working sets according to initial guesses.";
    case qpOASES::RET_MAX_NWSR_REACHED: return "Maximum number of working set recalculations performed.";
    case qpOASES::RET_CONSTRAINTS_NOT_SPECIFIED: return "Problem does comprise constraints! You also have to specify new constraints' bounds.";
    case qpOASES::RET_INVALID_FACTORISATION_FLAG: return "Invalid factorisation flag.";
    case qpOASES::RET_UNABLE_TO_SAVE_QPDATA: return "Unable to save QP data.";
    case qpOASES::RET_STEPDIRECTION_FAILED_TQ: return "Abnormal termination due to TQ factorisation.";
    case qpOASES::RET_STEPDIRECTION_FAILED_CHOLESKY: return "Abnormal termination due to Cholesky factorisation.";
    case qpOASES::RET_CYCLING_DETECTED: return "Cycling detected.";
    case qpOASES::RET_CYCLING_NOT_RESOLVED: return "Cycling cannot be resolved, QP probably infeasible.";
    case qpOASES::RET_CYCLING_RESOLVED: return "Cycling probably resolved.";
    case qpOASES::RET_STEPSIZE: return "For displaying performed stepsize.";
    case qpOASES::RET_STEPSIZE_NONPOSITIVE: return "For displaying non-positive stepsize.";
    case qpOASES::RET_SETUPSUBJECTTOTYPE_FAILED: return "Setup of SubjectToTypes failed.";
    case qpOASES::RET_ADDCONSTRAINT_FAILED: return "Addition of constraint to working set failed.";
    case qpOASES::RET_ADDCONSTRAINT_FAILED_INFEASIBILITY: return "Addition of constraint to working set failed (due to QP infeasibility).";
    case qpOASES::RET_ADDBOUND_FAILED: return "Addition of bound to working set failed.";
    case qpOASES::RET_ADDBOUND_FAILED_INFEASIBILITY: return "Addition of bound to working set failed (due to QP infeasibility).";
    case qpOASES::RET_REMOVECONSTRAINT_FAILED: return "Removal of constraint from working set failed.";
    case qpOASES::RET_REMOVEBOUND_FAILED: return "Removal of bound from working set failed.";
    case qpOASES::RET_REMOVE_FROM_ACTIVESET: return "Removing from active set...";
    case qpOASES::RET_ADD_TO_ACTIVESET: return "Adding to active set...";
    case qpOASES::RET_REMOVE_FROM_ACTIVESET_FAILED: return "Removing from active set failed.";
    case qpOASES::RET_ADD_TO_ACTIVESET_FAILED: return "Adding to active set failed.";
    case qpOASES::RET_CONSTRAINT_ALREADY_ACTIVE: return "Constraint is already active.";
    case qpOASES::RET_ALL_CONSTRAINTS_ACTIVE: return "All constraints are active, no further constraint can be added.";
    case qpOASES::RET_LINEARLY_DEPENDENT: return "New bound/constraint is linearly dependent.";
    case qpOASES::RET_LINEARLY_INDEPENDENT: return "New bound/constraint is linearly independent.";
    case qpOASES::RET_LI_RESOLVED: return "Linear indepence of active contraint matrix successfully resolved.";
    case qpOASES::RET_ENSURELI_FAILED: return "Failed to ensure linear indepence of active contraint matrix.";
    case qpOASES::RET_ENSURELI_FAILED_TQ: return "Abnormal termination due to TQ factorisation.";
    case qpOASES::RET_ENSURELI_FAILED_NOINDEX: return "QP is infeasible.";
    case qpOASES::RET_ENSURELI_FAILED_CYCLING: return "QP is infeasible.";
    case qpOASES::RET_BOUND_ALREADY_ACTIVE: return "Bound is already active.";
    case qpOASES::RET_ALL_BOUNDS_ACTIVE: return "All bounds are active, no further bound can be added.";
    case qpOASES::RET_CONSTRAINT_NOT_ACTIVE: return "Constraint is not active.";
    case qpOASES::RET_BOUND_NOT_ACTIVE: return "Bound is not active.";
    case qpOASES::RET_HESSIAN_NOT_SPD: return "Projected Hessian matrix not positive definite.";
    case qpOASES::RET_HESSIAN_INDEFINITE: return "Hessian matrix is indefinite.";
    case qpOASES::RET_MATRIX_SHIFT_FAILED: return "Unable to update matrices or to transform vectors.";
    case qpOASES::RET_MATRIX_FACTORISATION_FAILED: return "Unable to calculate new matrix factorisations.";
    case qpOASES::RET_PRINT_ITERATION_FAILED: return "Unable to print information on current iteration.";
    case qpOASES::RET_NO_GLOBAL_MESSAGE_OUTPUTFILE: return "No global message output file initialised.";
    case qpOASES::RET_DISABLECONSTRAINTS_FAILED: return "Unable to disbable constraints.";
    case qpOASES::RET_ENABLECONSTRAINTS_FAILED: return "Unable to enbable constraints.";
    case qpOASES::RET_ALREADY_ENABLED: return "Bound or constraint is already enabled.";
    case qpOASES::RET_ALREADY_DISABLED: return "Bound or constraint is already disabled.";
    case qpOASES::RET_NO_HESSIAN_SPECIFIED: return "No Hessian matrix has been specified.";
    case qpOASES::RET_USING_REGULARISATION: return "Using regularisation as Hessian matrix is not positive definite.";
    case qpOASES::RET_EPS_MUST_BE_POSITVE: return "Eps for regularisation must be sufficiently positive.";
    case qpOASES::RET_REGSTEPS_MUST_BE_POSITVE: return "Maximum number of regularisation steps must be non-negative.";
    case qpOASES::RET_HESSIAN_ALREADY_REGULARISED: return "Hessian has been already regularised.";
    case qpOASES::RET_CANNOT_REGULARISE_IDENTITY: return "Identity Hessian matrix cannot be regularised.";
    case qpOASES::RET_NO_REGSTEP_NWSR: return "No additional regularisation step could be performed due to limits.";
    case qpOASES::RET_FEWER_REGSTEPS_NWSR: return "Fewer additional regularisation steps have been performed due to limits.";
    case qpOASES::RET_CHOLESKY_OF_ZERO_HESSIAN: return "Cholesky decomposition of (unregularised) zero Hessian matrix.";
    case qpOASES::RET_CONSTRAINTS_ARE_NOT_SCALED: return "When defining __MANY_CONSTRAINTS__, l1 norm of each constraint must be not greater than one.";
    case qpOASES::RET_ERROR_IN_CONSTRAINTPRODUCT: return "Error in user-defined constraint product function.";
    // SQProblem
    case qpOASES::RET_UPDATEMATRICES_FAILED: return "Unable to update QP matrices.";
    case qpOASES::RET_UPDATEMATRICES_FAILED_AS_QP_NOT_SOLVED: return "Unable to update matrices as previous QP is not solved.";
    // Utils
    case qpOASES::RET_UNABLE_TO_OPEN_FILE: return "Unable to open file.";
    case qpOASES::RET_UNABLE_TO_WRITE_FILE: return "Unable to write into file.";
    case qpOASES::RET_UNABLE_TO_READ_FILE: return "Unable to read from file.";
    case qpOASES::RET_FILEDATA_INCONSISTENT: return "File contains inconsistent data.";
    // SolutionAnalysis
    case qpOASES::RET_UNABLE_TO_ANALYSE_QPROBLEM: return "Unable to analyse (S)QProblem(B) object";
    // Benchmark
    case qpOASES::RET_NWSR_SET_TO_ONE: return "Maximum number of working set changes was set to 1.";
    case qpOASES::RET_BENCHMARK_ABORTED: return "Benchmark aborted.";
    case qpOASES::RET_UNABLE_TO_READ_BENCHMARK: return "Unable to read benchmark data.";
    case qpOASES::RET_INITIAL_QP_SOLVED: return "Initial QP solved.";
    case qpOASES::RET_QP_SOLUTION_STARTED: return "Solving QP...";
    case qpOASES::RET_BENCHMARK_SUCCESSFUL: return "Benchmark terminated successfully.";
  }
  
  // Default error message
  stringstream ss;
  ss << "Unknown error flag: " << flag << ". Consult qpOASES documentation.";
  return ss.str();
}


} // namespace Interfaces
} // namespace CasADi

