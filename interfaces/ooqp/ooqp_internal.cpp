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

#include "ooqp_internal.hpp"

#include "casadi/matrix/sparsity_tools.hpp"
#include "casadi/matrix/matrix_tools.hpp"
#include "casadi/stl_vector_tools.hpp"

#include "casadi/pre_c99_support.hpp"

#include "Status.h"
// #define getPtr(x) (&(x).front())

using namespace std;
namespace CasADi {
namespace Interfaces {

OOQPInternal* OOQPInternal::clone() const{
  // Return a deep copy
  OOQPInternal* node = new OOQPInternal(input(QP_H).sparsity(),input(QP_A).sparsity());
  if(!node->is_init_)
    node->init();
  return node;
}
  
OOQPInternal::OOQPInternal(const CRSSparsity& H, const CRSSparsity& A) : QPSolverInternal(H,A){
  addOption("print_level",OT_INTEGER,0,"Print level. OOQP listends to print_level 0, 10 and 100");
  addOption("mutol",OT_REAL,1e-8,"tolerance as provided with setMuTol to OOQP");
  addOption("artol",OT_REAL,1e-8,"tolerance as provided with setArTol to OOQP");
  
  qp_=0;
  prob_=0;
  vars_=0; 
  resid_=0;
  s_=0;
}

OOQPInternal::~OOQPInternal(){ 

  if (s_) delete s_;
  if (qp_) {
    delete qp_;
    delete prob_;
    delete vars_;
    delete resid_;
  }
}

void OOQPInternal::evaluate(int nfdir, int nadir) {
  casadi_assert_message(nfdir==0 && nadir==0, "OOQPSolve::evaluate() not implemented for forward or backward mode");
  if (qp_==0) {
    allocate();
    casadi_assert(qp_!=0);
  }
  
  // Access nonzeros
  const vector<double> &data_A = input(QP_A).data();
  vector<double> &data_Aeq = A_eq_.data();
  vector<double> &data_Aineq = A_ineq_.data();
  
  // Get bounds 
  const vector<double>& lba = input(QP_LBA).data();
  const vector<double>& uba = input(QP_UBA).data();
  
  // Sparsity of A
  const vector<int> &rowind_A = input(QP_A).sparsity().rowind();

  // Loop over the rows of A corresponding to equalities
  const vector<int> &rowind_Aeq = A_eq_.sparsity().rowind();
  for(int i=0; i<eq_.size(); ++i){
    // Row of A
    int iA = eq_[i];
    
    // Copy the nonzeros of A
    int el, elA;
    for(el=rowind_Aeq[i], elA=rowind_A[iA]; el<rowind_Aeq[i+1]; ++el, ++elA){
      data_Aeq[el] = data_A[elA];
    }
    
    // Copy bounds
    bA_eq_[i] = lba[iA];
  }

  // Loop over the rows of A corresponding to inequalities
  const vector<int> &rowind_Aineq = A_ineq_.sparsity().rowind();
  for(int i=0; i<ineq_.size(); ++i){
    // Row of A
    int iA = ineq_[i];
    
    // Copy the nonzeros of A
    int el, elA;
    for(el=rowind_Aineq[i], elA=rowind_A[iA]; el<rowind_Aineq[i+1]; ++el, ++elA){
      data_Aineq[el] = data_A[elA];
    }
    
    // Copy bounds
    lbA_ineq_[i] = iclow_[i] ? lba[iA] : 0;
    ubA_ineq_[i] = icupp_[i] ? uba[iA] : 0;
  }
  
  // Lower variable bounds, 0 if infinity
  input(QP_LBX).get(lbX_);
  for(int k=0; k<lbX_.size();++k)
    if(ixlow_[k]==0)
      lbX_[k] = 0;

  // Upper variable bounds, 0 if infinity
  input(QP_UBX).get(ubX_);
  for(int k=0; k<ubX_.size();++k)
    if(ixupp_[k]==0)
      ubX_[k] = 0;
  
  // Hessian data
  input(QP_H).get(Hl_.data(),SPARSESYM);
  
  // Just calling s_->solve repeatedly on a GondzioSolver is a memory leak
  // So we must always allocate a fresh solver
  if (s_) delete s_;
  s_ = new GondzioSolver( qp_, prob_ );
  
  s_->setMuTol(double(getOption("mutol")));
  s_->setArTol(double(getOption("artol")));
  
  int flag = s_->solve(prob_,vars_, resid_);
    
  // Get primal solution
  casadi_assert(vars_->x->length()==output(QP_PRIMAL).size());
  vars_->x->copyIntoArray(getPtr(output(QP_PRIMAL)));
  if (isnan(output(QP_PRIMAL).at(0))) {
    casadi_warning("nan in decision variables. You probably need to do a solver.reInit() call before evaluate().");
  }
  
  // Get multipliers for the bounds
  vector<double> &dual_x = output(QP_DUAL_X).data();

  // BUG?
  #if 0
  // Lower bounds
  casadi_assert(dual_x.size()==vars_->gamma->length());
  vars_->gamma->copyIntoArray(getPtr(dual_x));
  
  // Upper bounds
  casadi_assert(dual_x.size()==vars_->phi->length());
  vars_->phi->copyIntoArray(getPtr(temp_));
    for(int k=0; k<dual_x.size(); ++k){
      dual_x[k] -= temp_[k];
    }
  #else
  // Workaround not counting the infinite values?
  
  // Lower bounds
  int k_gamma=0;
  vars_->gamma->copyIntoArray(getPtr(temp_));
  for(int k=0; k<dual_x.size(); ++k){
    if(ixlow_[k]){
      dual_x[k] = temp_[k_gamma++];
    } else {
      dual_x[k] = 0;
    }
  }
  casadi_assert(k_gamma==vars_->gamma->length());
  
  // Upper bounds
  int k_phi=0;
  vars_->phi->copyIntoArray(getPtr(temp_));
  for(int k=0; k<dual_x.size(); ++k){
    if(ixupp_[k]){
      dual_x[k] -= temp_[k_phi++];
    }
  }
  casadi_assert(k_phi==vars_->phi->length());
  #endif

  
  // Get multipliers for the equality constraints
  vector<double> &dual_a = output(QP_DUAL_A).data();
  casadi_assert(vars_->y->length()==eq_.size());
  vars_->y->copyIntoArray(getPtr(temp_));
  for(int k=0; k<eq_.size(); ++k){
    dual_a[eq_[k]] = temp_[k];
  }
  
  // Get multipliers for the inequality constraints
  casadi_assert(vars_->z->length()==ineq_.size());
  vars_->z->copyIntoArray(getPtr(temp_));
  for(int k=0; k<ineq_.size(); ++k){
    dual_a[ineq_[k]] = temp_[k];
  }
  
  if(flag!=SUCCESSFUL_TERMINATION) ooqp_error("Solve",flag);

  
  output(QP_COST).at(0) = prob_->objectiveValue(vars_);
}

void OOQPInternal::allocate() {
  if (qp_) {
    delete qp_;
    delete prob_;
    delete vars_;
    delete resid_;
  }

  // Decide if constraints are equality or inequality
  ineq_.clear();
  eq_.clear();
  const vector<double>& uba = input(QP_UBA).data();
  const vector<double>& lba = input(QP_LBA).data();
  for(int k=0; k<uba.size(); ++k){
    if(uba[k]==lba[k]){
      eq_.push_back(k);
    } else {
      ineq_.push_back(k);
    }
  }
  
  // Set up a Sparse solver
  qp_ = new QpGenSparseMa27( nx_, eq_.size(), ineq_.size(), input(QP_H).size() , input(QP_G).size(), input(QP_A).size() );
  

  // Split A in equalities and inequalities
  vector<int> all_A = range(0,input(QP_A).size2());
  A_eq_   = input(QP_A)(eq_,all_A);
  A_ineq_ = input(QP_A)(ineq_,all_A);
  
  // Split up LBA & UBA in equalities and inequalities
  bA_eq_.resize(eq_.size());
  
  lbA_ineq_ = input(QP_LBA)(ineq_,0).data();
  ubA_ineq_ = input(QP_UBA)(ineq_,0).data();
  
  // Copy the decision variables bounds
  lbX_.resize(nx_);
  ubX_.resize(nx_);
  
  // Infinities on LBX & UBX are set to zero (required for OOQp)
  for (int k=0; k<input(QP_LBX).size(); ++k) {
    if (input(QP_LBX).at(k)==-numeric_limits<double>::infinity()) 
      lbX_.at(k)=0;
  }
  for (int k=0; k<input(QP_UBX).size(); ++k) {
    if (input(QP_UBX).at(k)==numeric_limits<double>::infinity()) 
      ubX_.at(k)=0;
  }
  
  
  // Set ixlow & ixupp: they have to be zero for infinite bounds and 1 otherwise
  for (int k=0; k<input(QP_LBX).size(); ++k) 
    ixlow_[k] = (input(QP_LBX).at(k) != -numeric_limits<double>::infinity());
  for (int k=0; k<input(QP_UBX).size(); ++k) 
    ixupp_[k] = (input(QP_UBX).at(k) != numeric_limits<double>::infinity());
  
  // Set iclow & icupp: they have to be zero for infinite bounds and 1 otherwise
  iclow_.resize(ineq_.size(),0);
  icupp_.resize(ineq_.size(),0);
  
  for (int k=0; k<lbA_ineq_.size();++k) 
    iclow_[k] = (lbA_ineq_.at(k) != -numeric_limits<double>::infinity());
  
  for (int k=0; k<ubA_ineq_.size();++k) 
    icupp_[k] = (ubA_ineq_.at(k) != numeric_limits<double>::infinity());

  // Infinities on LBA & UBA are set to zero (required for OOQp)
  for (int k=0; k<lbA_ineq_.size();k++) {
    if (lbA_ineq_.at(k) == -numeric_limits<double>::infinity()) 
      lbA_ineq_.at(k)=0;
  } 
  
  for (int k=0; k<ubA_ineq_.size();k++) {
    if (ubA_ineq_.at(k) == numeric_limits<double>::infinity()) 
      ubA_ineq_.at(k)=0;
  }
  
  // Get references to sparsity (NOTE: OOQP does not appear do be const correct)
  int *Hl_rowind = const_cast<int*>(getPtr(Hl_.sparsity().rowind()));
  int *Hl_col = const_cast<int*>(getPtr(Hl_.sparsity().col()));

  int *eq_rowind = const_cast<int*>(getPtr(A_eq_.sparsity().rowind()));
  int *eq_col = const_cast<int*>(getPtr(A_eq_.sparsity().col()));
  
  int *ineq_rowind = const_cast<int*>(getPtr(A_ineq_.sparsity().rowind()));
  int *ineq_col = const_cast<int*>(getPtr(A_ineq_.sparsity().col()));
  
  // Set all pointers to problem data
  prob_ = (QpGenData * )qp_->makeData( getPtr(input(QP_G)),
                                       Hl_rowind,  Hl_col,  getPtr(Hl_),
                                       getPtr(lbX_),  getPtr(ixlow_),
                                       getPtr(ubX_),  getPtr(ixupp_),
                                       eq_rowind, eq_col,  getPtr(A_eq_),
                                       getPtr(bA_eq_),
                                       ineq_rowind, ineq_col,  getPtr(A_ineq_),
                                       getPtr(lbA_ineq_),  getPtr(iclow_),
                                       getPtr(ubA_ineq_),  getPtr(icupp_));

  // Further setup of the QP problem
  vars_ = (QpGenVars *) qp_->makeVariables( prob_ );
  resid_ = (QpGenResiduals *) qp_->makeResiduals( prob_ );
}

void OOQPInternal::init(){
   // Call the init method of the base class
  QPSolverInternal::init();
  
  // Structure to hold the sparsity and data of the lower triangular part of the Hessian
  Hl_ = DMatrix(lowerSparsity(input(QP_H).sparsity()));
  
  // Allocate vectors to hold indicators of infinite variable bounds
  ixlow_.resize(nx_,0);
  ixupp_.resize(nx_,0);

  // Temporary vector
  temp_.resize(std::max(nx_,nc_));
}

map<int,string> OOQPInternal::calc_flagmap(){
  map<int,string> f;

  f[SUCCESSFUL_TERMINATION] = "SUCCESSFUL_TERMINATION";
  f[NOT_FINISHED] = "NOT_FINISHED";
  f[MAX_ITS_EXCEEDED] = "MAX_ITS_EXCEEDED";
  f[INFEASIBLE] = "INFEASIBLE";
  f[UNKNOWN] = "UNKNOWN";
  return f;
}
  
map<int,string> OOQPInternal::flagmap = OOQPInternal::calc_flagmap();

void OOQPInternal::ooqp_error(const string& module, int flag){
  // Find the error
  map<int,string>::const_iterator it = flagmap.find(flag);
  
  stringstream ss;
  if(it == flagmap.end()){
    ss << "Unknown error (" << flag << ") from module \"" << module << "\".";
  } else {
    ss << "Module \"" << module << "\" returned flag \"" << it->second << "\".";
  }
  ss << " Consult OOQP documentation.";
  casadi_warning(ss.str());
}

} // namespace Interfaces
} // namespace CasADi

// #undef getPtr
