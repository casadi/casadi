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

CplexInternal::CplexInternal(const CRSSparsity& H, const CRSSparsity& A) : QPSolverInternal(H, A){
  // Options available
  addOption("qp_method",    OT_INTEGER, 0, "Determines which CPLEX algorithm to use. \
    0: Automatic, \
    1: Primal simplex, \
    2: Dual simplex, \
    3: Network, \
    4: Barrier, \
    5: Sifting, \
    6: Concurent. \
    7: Crossover (start with barrier and use simplex later on with warm-start)\
    Default: 0");
  addOption("dump_to_file",   OT_BOOLEAN,        false, "Dumps QP to file in CPLEX format. Default: false");
  addOption("dump_filename",   OT_STRING,     "qp.dat", "The filename to dump to. Default: qp.dat");
  addOption("debug",          OT_BOOLEAN,        false, "Print debug information");
  addOption("tol",               OT_REAL,         1E-6, "Tolerance of solver");
  addOption("dep_check",      OT_INTEGER,            0, "Detect redundant constraints. \
    -1: automatic\
     0: off\
     1: at the beginning of preprocessing\
     2: at the end of preprocessing\
     3: at the begining and at the end of preprocessing");
  addOption("simplex_maxiter", OT_INTEGER,       2100000000, "Maximum number of simplex iterations.");
  addOption("barrier_maxiter", OT_INTEGER,       2100000000, "Maximum number of barrier iterations.");
  addOption("warm_start",      OT_BOOLEAN,            false, "Use warm start with simplex methods (affects only the simplex methods).");
  addOption("convex",          OT_BOOLEAN,             true, "Indicates if the QP is convex or not (affects only the barrier method).");
  
  // Initializing members
  // Number of vars
  NUMCOLS_ = H.size1();
  // Number of constraints
  NUMROWS_ = A.size1();
  // Setting warm-start flag
  is_warm_ = false;
}
void CplexInternal::init(){
  qp_method_     = getOption("qp_method");
  if (qp_method_ < 0 || qp_method_ > 7){
    casadi_error("Invalid QP method given.");
  }
  dump_to_file_  = getOption("dump_to_file");
  debug_ = getOption("debug");
  tol_ = getOption("tol");
//  dump_filename_ = getOption("dump_filename");
  debug_ = getOption("debug");

  int status;
  env_ = 0;
  QPSolverInternal::init();
  env_ = CPXopenCPLEX (&status);
  if (!env_){
    cout << "CPLEX: Cannot initialize CPLEX environment. STATUS: " << status << "\n";
    throw CasadiException("Cannot initialize CPLEX environment.\n");
  }
  // Turn on some debug messages if requested
  if (debug_){
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
  status = CPXsetintparam(env_, CPX_PARAM_DEPIND, getOption("dep_check"));
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

  // Matrix A
  matbeg_.resize(NUMCOLS_);
  matcnt_.resize(NUMCOLS_);
  matind_.resize(input(QP_A).size());
  matval_.resize(input(QP_A).size());
  // Matrix H
  qmatbeg_.resize(NUMCOLS_);
  qmatcnt_.resize(NUMCOLS_);
  qmatind_.resize(input(QP_H).size());
  qmatval_.resize(input(QP_H).size());
  
  // Linear term
  if (!isDense(input(QP_G))){
    casadi_error("input(QP_G) must be dense.");
  }
  if (!isDense(input(QP_LBA))){
    casadi_error("input(QP_LBA) must be dense.");
  }
  if (!isDense(input(QP_UBA))){
    casadi_error("input(QP_UBA) must be dense.");
  }
  if (!isDense(input(QP_LBX))){
    casadi_error("input(QP_LBX) must be dense.");
  }
  if (!isDense(input(QP_UBX))){
    casadi_error("input(QP_UBX) must be dense.");
  }

  lp_ = CPXcreateprob(env_, &status, "QP from CasADi");
}

void CplexInternal::evaluate(int nfdir, int nadir){

  casadi_assert(nfdir == 0 && nadir == 0);

  int status;

  // We change method in crossover
  if ( is_warm_ && qp_method_ == 7){
    status = CPXsetintparam(env_, CPX_PARAM_QPMETHOD, 1);
  }

  obj_ = input(QP_G).ptr();
  lb_ = input(QP_LBX).ptr();
  ub_ = input(QP_UBX).ptr();

  // Looping over constraints
  for(int i = 0; i < NUMROWS_; ++i){
    // CPX_INFBOUND
  
    // Equality
    if (input(QP_UBA).elem(i) - input(QP_LBA).elem(i) < 1E-20){
      sense_[i] = 'E';
      rhs_[i] = input(QP_LBA).elem(i);
      rngval_[i] = 0.;
    }
    // Ineq - no lower bound
    else if (input(QP_LBA).elem(i) < -CPX_INFBOUND){
      sense_[i] = 'L';
      rhs_[i] = input(QP_UBA).elem(i);
      //rngval_[i] = input(QP_UBA).elem(i) - input(QP_LBA).elem(i);
      rngval_[i] = 0.;
    }
    // Ineq - no upper bound
    else if (input(QP_UBA).elem(i) > CPX_INFBOUND){
      sense_[i] = 'G';
      rhs_[i] = input(QP_LBA).elem(i);
      rngval_[i] = 0.;
    }
    // Inew both upper and lower bounds
    else{
      sense_[i] = 'R';
      rhs_[i] = input(QP_LBA).elem(i);
      rngval_[i] = input(QP_UBA).elem(i) - input(QP_LBA).elem(i);
    }
  }

  // Preparing coefficient matrix A
  dmatrixToCplex(input(QP_A), matbeg_.data(), matcnt_.data(), matind_.data(), matval_.data());
  // Copying objective, constraints, and bounds.
  status = CPXcopylp (env_, lp_, NUMCOLS_, NUMROWS_, objsen_, obj_, rhs_.data(),
       sense_.data(), matbeg_.data(), matcnt_.data(), matind_.data(), matval_.data(), lb_, ub_, rngval_.data()); 

  // Preparing coefficient matrix Q
  dmatrixToCplex(input(QP_H), qmatbeg_.data(), qmatcnt_.data(), qmatind_.data(), qmatval_.data());
  // Copying quadratic term
  status = CPXcopyquad(env_, lp_, qmatbeg_.data(), qmatcnt_.data(), qmatind_.data(), qmatval_.data());

  if (dump_to_file_){
    const char* fn = string(getOption("dump_filename")).c_str();
    CPXwriteprob(env_, lp_, fn, "LP");
  }

  // Warm-starting if possible
  if (qp_method_ != 0 && qp_method_ != 4 && is_warm_){
    // TODO: Initialize slacks and dual variables of bound constraints
    CPXcopystart(env_, lp_, cstat_.data(), rstat_.data(), input(QP_X_INIT).ptr(), NULL, NULL, input(QP_LAMBDA_INIT).ptr());
  }
  else{
    status = CPXcopystart(env_, lp_, NULL, NULL, input(QP_X_INIT).ptr(), NULL, NULL, input(QP_LAMBDA_INIT).ptr());
  }

  // Optimize...
  status = CPXqpopt(env_, lp_);

  if (status){
    casadi_error("CPLEX: Failed to solve QP...");
  }
  // Retrieving solution
  int solstat; 
  int objval;
  std::vector<double> slack;
  slack.resize(NUMROWS_);
  status = CPXsolution (env_, lp_, &solstat,
   output(QP_COST).ptr(), 
   output(QP_PRIMAL).ptr(), 
   output(QP_LAMBDA_A).ptr(),
   slack.data(),
   output(QP_LAMBDA_X).ptr()
  ); 
  
  if(status){
    cout << "CPLEX: Failed to get solution.\n";
  } 
  // Retrieving the basis
  if (qp_method_ != 0 && qp_method_ != 4){
    status = CPXgetbase(env_, lp_, cstat_.data(), rstat_.data());
  }

  // Flip the sign of the multipliers
  for (int k=0;k<output(QP_LAMBDA_A).size();++k) output(QP_LAMBDA_A).data()[k]= - output(QP_LAMBDA_A).data()[k];
  for (int k=0;k<output(QP_LAMBDA_X).size();++k) output(QP_LAMBDA_X).data()[k]= - output(QP_LAMBDA_X).data()[k];
  
  int solnstat = CPXgetstat (env_, lp_);
  string errormsg;
  if(debug_){
    if      (solnstat == CPX_STAT_OPTIMAL){
      errormsg = string("CPLEX: solution status: Optimal solution found.\n");
    }
    else if (solnstat == CPX_STAT_UNBOUNDED) {
      errormsg = string("CPLEX: solution status: Model is unbounded\n");
    }
    else if (solnstat == CPX_STAT_INFEASIBLE) {
      errormsg = string("CPLEX: solution status: Model is infeasible\n");
    }
    else if (solnstat == CPX_STAT_INForUNBD) {
      errormsg = string("CPLEX: solution status: Model is infeasible or unbounded\n");
    }
    else if (solnstat == CPX_STAT_OPTIMAL_INFEAS){
      errormsg = string("CPLEX: solution status: Optimal solution is available but with infeasibilities\n");
    }
    else if (solnstat == CPX_STAT_NUM_BEST){
      errormsg = string("CPLEX: solution status: Solution available, but not proved optimal due to numeric difficulties.\n");
    }
    else if (solnstat == CPX_STAT_FIRSTORDER){
      errormsg = string("CPLEX: solution status: Solution satisfies first-order optimality conditions, but is not necessarily globally optimal.\n");
    }
    else{
      errormsg = string("CPLEX: solution status: ") + to_string(solnstat) + string("\n");
    }
    cout << errormsg;

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
  CplexInternal* node = new CplexInternal(input(QP_H).sparsity(), input(QP_A).sparsity());
  if(!node->is_init_)
    node->init();
  return node;
}

/// Destructor
CplexInternal::~CplexInternal(){
  int status; 
  status = CPXfreeprob (env_, &lp_);
  if ( status ) {
     std::cout << "CPXfreeprob failed, error code " << status << ".\n";
  }
  // Closing down license
  status = CPXcloseCPLEX (&env_);
}
void CplexInternal::dmatrixToCplex(DMatrix& M, int* matbeg, int* matcnt, int* matind, double* matval){
  std::vector<int> rowind,col;

  DMatrix MT = M.trans();
//  DMatrix MT = M;
  MT.sparsity().getSparsityCRS(rowind, col);
  std::vector<double>& data = MT.data();
//  int i, j;
//  int j = 0;
//  double val;
  int nnz_counter = 0;
  int nnz_col = 0;
  for(int r=0; r<rowind.size()-1; ++r){
    matbeg[r] = nnz_counter;
    for(int el=rowind[r]; el<rowind[r+1]; ++el){
//      i = r;
//      j = col[el];
//      val = data[el];
      matind[nnz_counter] = col[el];
      matval[nnz_counter] = data[el];
      ++nnz_counter;
      ++nnz_col;
    }
    matcnt[r] = nnz_col;
    nnz_col = 0;
  }
};
} // end namespace CasADi

