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
    Default: 0");
  addOption("dump_to_file",   OT_BOOLEAN,        false, "Dumps QP to file in CPLEX format. Default: false");
  addOption("dump_filename",   OT_STRING, "qp.dat", "The filename to dump to. Default: qp.dat");
  addOption("debug",          OT_BOOLEAN,        false, "Print debug information");
  addOption("tol",             OT_REAL,     1E-6, "Tolerance of solver");
  
  // Initializing members
  // Number of vars
  NUMCOLS_ = H.size1();
  // Number of constraints
  NUMROWS_ = A.size1();
}
void CplexInternal::init(){
  qp_method_     = getOption("qp_method");
  dump_to_file_  = getOption("dump_to_file");
  debug_ = getOption("debug");
  tol_ = getOption("tol");
//  dump_filename_ = getOption("dump_filename");
  debug_ = getOption("debug");

  int status;
  env_ = 0;
  QPSolverInternal::init();
//  CPXsetintparam (env_, CPX_PARAM_SCRIND, CPX_OFF);
  env_ = CPXopenCPLEX (&status);
  if (!env_){
    std::cerr << "Cannot initialize CPLEX environment. STATUS: " << status << "\n";
    throw("Cannot initialize CPLEX environment.\n");
  }
  // Turn on some debug messages if requested
  if (debug_){
    CPXsetintparam (env_, CPX_PARAM_SCRIND, CPX_ON);
  }
  else{
    CPXsetintparam (env_, CPX_PARAM_SCRIND, CPX_OFF);
  }
  if (status){
    std::cout << "Problem with setting parameter... ERROR: " << status << std::endl;
  }


  // Optimality tolerance
  CPXsetdblparam(env_, CPX_PARAM_EPOPT, tol_);
  // Feasibility tolerance
  CPXsetdblparam(env_, CPX_PARAM_EPRHS, tol_);
  
  // Doing allocation of CPLEX data
  // Objective is to be minimized
  objsen_ = CPX_MIN;
  // Allocation of data
  sense_.resize(NUMROWS_);
  rhs_.resize(NUMROWS_);
  rngval_.resize(NUMROWS_);

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

}

void CplexInternal::evaluate(int nfdir, int nadir){

  int status;

  casadi_assert(nfdir == 0 && nadir == 0);

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

    lp_ = CPXcreateprob(env_, &status, "QP from CasADi");
  }
  
  // Preparing coefficient matrix A
  dmatrixToCplex(input(QP_A), matbeg_.data(), matcnt_.data(), matind_.data(), matval_.data());
  // Copying objective, constraints, and bounds.
  status = CPXcopylp (env_, lp_, NUMCOLS_, NUMROWS_, objsen_, obj_, rhs_.data(),
       sense_.data(), matbeg_.data(), matcnt_.data(), matind_.data(), matval_.data(), lb_, ub_, rngval_.data()); 

  // Preparing coefficient matrix Q
  dmatrixToCplex(input(QP_H), qmatbeg_.data(), qmatcnt_.data(), qmatind_.data(), qmatval_.data());

  status = CPXcopyquad(env_, lp_, qmatbeg_.data(), qmatcnt_.data(), qmatind_.data(), qmatval_.data());

  // Setting optimization method
  status = CPXsetintparam(env_, CPX_PARAM_QPMETHOD, qp_method_);
  

  if (dump_to_file_){
    CPXwriteprob(env_, lp_, "qp.lp", NULL);
  }

  // Optimize...
  status = CPXqpopt(env_, lp_);

  if (status){
    casadi_error("Failed to solve QP...");
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
  
  // Flip the sign of the multipliers
  for (int k=0;k<output(QP_LAMBDA_A).size();++k) output(QP_LAMBDA_A).data()[k]= - output(QP_LAMBDA_A).data()[k];
  for (int k=0;k<output(QP_LAMBDA_X).size();++k) output(QP_LAMBDA_X).data()[k]= - output(QP_LAMBDA_X).data()[k];
  
  if(status){
    std::cerr << "CPLEX: Failed to get solution.\n";
  } 
  if(debug_){
    int solnstat = CPXgetstat (env_, lp_);
    if      ( solnstat == CPX_STAT_UNBOUNDED ) {
      printf ("Model is unbounded\n");
    }
    else if ( solnstat == CPX_STAT_INFEASIBLE ) {
      printf ("Model is infeasible\n");
    }
    else if ( solnstat == CPX_STAT_INForUNBD ) {
      printf ("Model is infeasible or unbounded\n");
    }
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

