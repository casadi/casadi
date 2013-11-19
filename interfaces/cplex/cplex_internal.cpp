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
#include "cplex_internal.hpp"
#include "symbolic/stl_vector_tools.hpp"
#include "symbolic/matrix/matrix_tools.hpp"
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <vector>

#include "ilcplex/cplex.h"

namespace CasADi{

using namespace std;

CplexInternal::CplexInternal(const std::vector<CRSSparsity>& st) : QPSolverInternal(st){
  // Options available
  addOption("qp_method",    OT_STRING, "automatic", "Determines which CPLEX algorithm to use.","automatic|primal_simplex|dual_simplex|network|barrier|sifting|concurrent|crossover");
  addOption("dump_to_file",   OT_BOOLEAN,        false, "Dumps QP to file in CPLEX format.");
  addOption("dump_filename",   OT_STRING,     "qp.dat", "The filename to dump to.");
  addOption("tol",               OT_REAL,         1E-6, "Tolerance of solver");
  addOption("dep_check",      OT_STRING,         "off", "Detect redundant constraints.","automatic:-1|off:0|begin:1|end:2|both:3");
  addOption("simplex_maxiter", OT_INTEGER,       2100000000, "Maximum number of simplex iterations.");
  addOption("barrier_maxiter", OT_INTEGER,       2100000000, "Maximum number of barrier iterations.");
  addOption("warm_start",      OT_BOOLEAN,            false, "Use warm start with simplex methods (affects only the simplex methods).");
  addOption("convex",          OT_BOOLEAN,             true, "Indicates if the QP is convex or not (affects only the barrier method).");
  
  const CRSSparsity& A = st_[QP_STRUCT_A];
  const CRSSparsity& H = st_[QP_STRUCT_H];
  
  // Initializing members
  // Number of vars
  NUMCOLS_ = H.size1();
  // Number of constraints
  NUMROWS_ = A.size1();
  // Setting warm-start flag
  is_warm_ = false;

  // Set pointer to zero to avoid deleting a nonexisting instance
  env_ = 0;
  lp_ = 0;
}
void CplexInternal::init(){
  // Free any existing Cplex instance
  freeCplex();

  // Call the init method of the base class
  QPSolverInternal::init();
  
  qp_method_     = getOptionEnumValue("qp_method");
  dump_to_file_  = getOption("dump_to_file");
  tol_ = getOption("tol");
  //  dump_filename_ = getOption("dump_filename");
  
  int status;
  casadi_assert(env_==0);
  env_ = CPXopenCPLEX (&status);
  casadi_assert_message(env_!=0, "CPLEX: Cannot initialize CPLEX environment. STATUS: " << status);

  // Turn on some debug messages if requested
  if (verbose()){
    CPXsetintparam (env_, CPX_PARAM_SCRIND, CPX_ON);
  }
  else{
    CPXsetintparam (env_, CPX_PARAM_SCRIND, CPX_OFF);
  }
  if (status){
    std::cout << "CPLEX: Problem with setting parameter... ERROR: " << status << std::endl;
  }

  /* SETTING OPTIONS */
  // Optimality tolerance
  status = CPXsetdblparam(env_, CPX_PARAM_EPOPT, tol_);
  // Feasibility tolerance
  status = CPXsetdblparam(env_, CPX_PARAM_EPRHS, tol_);
  // We start with barrier if crossover was chosen.
  if (qp_method_ == 7){
    status = CPXsetintparam(env_, CPX_PARAM_QPMETHOD, 4);
    // Warm-start is default with this algorithm
    setOption("warm_start", true);
  }
  // Otherwise we just chose the algorithm
  else{
    status = CPXsetintparam(env_, CPX_PARAM_QPMETHOD, qp_method_);
  }
  // Setting dependency check option
  status = CPXsetintparam(env_, CPX_PARAM_DEPIND, getOptionEnumValue("dep_check"));
  // Setting barrier iteration limit
  status = CPXsetintparam(env_, CPX_PARAM_BARITLIM, getOption("barrier_maxiter"));
  // Setting simplex iteration limit
  status = CPXsetintparam(env_, CPX_PARAM_ITLIM, getOption("simplex_maxiter"));
  if (qp_method_ == 7){
    // Setting crossover algorithm
    status = CPXsetintparam(env_, CPX_PARAM_BARCROSSALG, 1);
  }
  if(!bool(getOption("convex"))){
    // Enabling non-convex QPs
    status = CPXsetintparam(env_, CPX_PARAM_SOLUTIONTARGET, CPX_SOLUTIONTARGET_FIRSTORDER);
  }

  // Exotic parameters, once they might become options...

  // Do careful numerics with numerically unstable problem
  //status = CPXsetintparam(env_, CPX_PARAM_NUMERICALEMPHASIS, 1);
  // Set scaling approach
  //status = CPXsetintparam(env_, CPX_PARAM_SCAIND, 1);
  // Set Markowitz tolerance
  //status = CPXsetdblparam(env_, CPX_PARAM_EPMRK, 0.9);
 
  // Doing allocation of CPLEX data
  // Objective is to be minimized
  objsen_ = CPX_MIN;

  // Allocation of data
  // Type of constraint
  sense_.resize(NUMROWS_);
  // Right-hand side of constraints
  rhs_.resize(NUMROWS_);
  // Range value for lower AND  upper bounded constraints
  rngval_.resize(NUMROWS_);
  // Basis for primal variables
  cstat_.resize(NUMCOLS_);
  rstat_.resize(NUMROWS_);

  // Matrix A (Cplex reqests its transpose)
  CRSSparsity AT_sparsity = input(QP_SOLVER_A).sparsity().transpose(AT_nonzero_mapping_);
  toCplexSparsity(AT_sparsity,matbeg_,matcnt_,matind_);
  matval_.resize(input(QP_SOLVER_A).size());
  
  // Matrix H
  toCplexSparsity(input(QP_SOLVER_H).sparsity(),qmatbeg_,qmatcnt_,qmatind_);
  
  // Linear term
  if (!isDense(input(QP_SOLVER_G))){
    casadi_error("input(QP_SOLVER_G) must be dense.");
  }
  if (!isDense(input(QP_SOLVER_LBA))){
    casadi_error("input(QP_SOLVER_LBA) must be dense.");
  }
  if (!isDense(input(QP_SOLVER_UBA))){
    casadi_error("input(QP_SOLVER_UBA) must be dense.");
  }
  if (!isDense(input(QP_SOLVER_LBX))){
    casadi_error("input(QP_SOLVER_LBX) must be dense.");
  }
  if (!isDense(input(QP_SOLVER_UBX))){
    casadi_error("input(QP_SOLVER_UBX) must be dense.");
  }
  
  casadi_assert(lp_==0);
  lp_ = CPXcreateprob(env_, &status, "QP from CasADi");
}

void CplexInternal::evaluate(){

  if (inputs_check_) checkInputs();

  int status;

  // We change method in crossover
  if ( is_warm_ && qp_method_ == 7){
    status = CPXsetintparam(env_, CPX_PARAM_QPMETHOD, 1);
  }

  obj_ = input(QP_SOLVER_G).ptr();
  lb_ = input(QP_SOLVER_LBX).ptr();
  ub_ = input(QP_SOLVER_UBX).ptr();

  // Looping over constraints
  for(int i = 0; i < NUMROWS_; ++i){
    // CPX_INFBOUND
  
    // Equality
    if (input(QP_SOLVER_UBA).elem(i) - input(QP_SOLVER_LBA).elem(i) < 1E-20){
      sense_[i] = 'E';
      rhs_[i] = input(QP_SOLVER_LBA).elem(i);
      rngval_[i] = 0.;
    }
    // Ineq - no lower bound
    else if (input(QP_SOLVER_LBA).elem(i) < -CPX_INFBOUND){
      sense_[i] = 'L';
      rhs_[i] = input(QP_SOLVER_UBA).elem(i);
      //rngval_[i] = input(QP_SOLVER_UBA).elem(i) - input(QP_SOLVER_LBA).elem(i);
      rngval_[i] = 0.;
    }
    // Ineq - no upper bound
    else if (input(QP_SOLVER_UBA).elem(i) > CPX_INFBOUND){
      sense_[i] = 'G';
      rhs_[i] = input(QP_SOLVER_LBA).elem(i);
      rngval_[i] = 0.;
    }
    // Inew both upper and lower bounds
    else{
      sense_[i] = 'R';
      rhs_[i] = input(QP_SOLVER_LBA).elem(i);
      rngval_[i] = input(QP_SOLVER_UBA).elem(i) - input(QP_SOLVER_LBA).elem(i);
    }
  }

  // Map the nonzeros of A to its transpose
  const vector<double>& A_data = input(QP_SOLVER_A).data();
  for(int k=0; k<A_data.size(); ++k){
    matval_[k] = A_data[AT_nonzero_mapping_[k]];
  }
  
  // Copying objective, constraints, and bounds.
  status = CPXcopylp (env_, lp_, NUMCOLS_, NUMROWS_, objsen_, obj_, rhs_.data(),
       sense_.data(), matbeg_.data(), matcnt_.data(), matind_.data(), matval_.data(), lb_, ub_, rngval_.data()); 

  // Preparing coefficient matrix Q
  const double* qmatval = input(QP_SOLVER_H).ptr();
  status = CPXcopyquad(env_, lp_, qmatbeg_.data(), qmatcnt_.data(), qmatind_.data(), qmatval);

  if (dump_to_file_){
    const char* fn = string(getOption("dump_filename")).c_str();
    CPXwriteprob(env_, lp_, fn, "LP");
  }

  // Warm-starting if possible
  if (qp_method_ != 0 && qp_method_ != 4 && is_warm_){
    // TODO: Initialize slacks and dual variables of bound constraints
    CPXcopystart(env_, lp_, cstat_.data(), rstat_.data(), input(QP_SOLVER_X0).ptr(), NULL, NULL, input(QP_SOLVER_LAM_X0).ptr());
  }
  else{
    status = CPXcopystart(env_, lp_, NULL, NULL, input(QP_SOLVER_X0).ptr(), NULL, NULL, input(QP_SOLVER_LAM_X0).ptr());
  }

  // Optimize...
  status = CPXqpopt(env_, lp_);

  if (status){
    casadi_error("CPLEX: Failed to solve QP...");
  }
  // Retrieving solution
  int solstat; 
  
  std::vector<double> slack;
  slack.resize(NUMROWS_);
  status = CPXsolution (env_, lp_, &solstat,
   output(QP_SOLVER_COST).ptr(), 
   output(QP_SOLVER_X).ptr(), 
   output(QP_SOLVER_LAM_A).ptr(),
   slack.data(),
   output(QP_SOLVER_LAM_X).ptr()
  ); 
  
  if(status){
    cout << "CPLEX: Failed to get solution.\n";
  } 
  // Retrieving the basis
  if (qp_method_ != 0 && qp_method_ != 4){
    status = CPXgetbase(env_, lp_, cstat_.data(), rstat_.data());
  }

  // Flip the sign of the multipliers
  for (int k=0;k<output(QP_SOLVER_LAM_A).size();++k) output(QP_SOLVER_LAM_A).data()[k]= - output(QP_SOLVER_LAM_A).data()[k];
  for (int k=0;k<output(QP_SOLVER_LAM_X).size();++k) output(QP_SOLVER_LAM_X).data()[k]= - output(QP_SOLVER_LAM_X).data()[k];
  
  int solnstat = CPXgetstat (env_, lp_);
  stringstream errormsg; // NOTE: Why not print directly to cout and cerr?
  if(verbose()){
    if      (solnstat == CPX_STAT_OPTIMAL){
      errormsg << "CPLEX: solution status: Optimal solution found.\n";
    }
    else if (solnstat == CPX_STAT_UNBOUNDED) {
      errormsg << "CPLEX: solution status: Model is unbounded\n";
    }
    else if (solnstat == CPX_STAT_INFEASIBLE) {
      errormsg << "CPLEX: solution status: Model is infeasible\n";
    }
    else if (solnstat == CPX_STAT_INForUNBD) {
      errormsg << "CPLEX: solution status: Model is infeasible or unbounded\n";
    }
    else if (solnstat == CPX_STAT_OPTIMAL_INFEAS){
      errormsg << "CPLEX: solution status: Optimal solution is available but with infeasibilities\n";
    }
    else if (solnstat == CPX_STAT_NUM_BEST){
      errormsg << "CPLEX: solution status: Solution available, but not proved optimal due to numeric difficulties.\n";
    }
    else if (solnstat == CPX_STAT_FIRSTORDER){
      errormsg << "CPLEX: solution status: Solution satisfies first-order optimality conditions, but is not necessarily globally optimal.\n";
    }
    else{
      errormsg << "CPLEX: solution status: " <<  solnstat << "\n";
    }
    cout << errormsg.str();

    // Printing basis condition number
    //double cn;
    //status = CPXgetdblquality(env_, lp_, &cn, CPX_KAPPA);
    //cout << "CPLEX: Basis condition number: " << cn << endl;
  }
  if (solnstat != CPX_STAT_OPTIMAL){
//    throw CasadiException(errormsg.c_str());
  }

  // Next time we warm start
  if (bool(getOption("warm_start"))){
    is_warm_ = true;
  }

}

CplexInternal* CplexInternal::clone() const{
  // Return a deepcopy
  CplexInternal* node = new CplexInternal(st_);
  if(!node->is_init_)
    node->init();
  return node;
}

CplexInternal::~CplexInternal(){
  freeCplex();
}

void CplexInternal::freeCplex(){
  // Return flag
  int status;

  // Only free if Cplex problem if it has been allocated
  if(lp_!=0){
    status = CPXfreeprob (env_, &lp_);
    if(status!=0){
      std::cerr << "CPXfreeprob failed, error code " << status << ".\n";
    }
    lp_ = 0;
  }
  
  // Closing down license
  if(env_!=0){
    status = CPXcloseCPLEX(&env_);
    if(status!=0){
      std::cerr << "CPXcloseCPLEX failed, error code " << status << ".\n";
    }        
    env_ = 0;
  }
}

void CplexInternal::toCplexSparsity(const CRSSparsity& sp_trans, vector<int> &matbeg, vector<int>& matcnt, vector<int>& matind){
  // Get sparsity
  int ncol = sp_trans.size1();
  //int nrow = sp_trans.size2();
  const std::vector<int>& colind = sp_trans.rowind();
  const std::vector<int>& row = sp_trans.col();

  // The row for each nonzero
  matind.resize(row.size());
  copy(row.begin(),row.end(),matind.begin());
  
  // The beginning of each column
  matbeg.resize(ncol);
  copy(colind.begin(),colind.begin()+ncol,matbeg.begin());
  
  // The number of elements in each column
  matcnt.resize(ncol);
  transform(colind.begin()+1,colind.end(),colind.begin(),matcnt.begin(),minus<int>());
}

} // end namespace CasADi

