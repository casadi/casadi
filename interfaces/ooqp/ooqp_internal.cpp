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

using namespace std;
namespace CasADi {
namespace Interfaces {

OOQPInternal* OOQPInternal::clone() const{
  // Return a deep copy
  OOQPInternal* node = new OOQPInternal(H,G,A);
  if(!node->is_init)
    node->init();
  return node;
}
  
OOQPInternal::OOQPInternal(const CRSSparsity & H, const CRSSparsity & G, const CRSSparsity & A) : QPSolverInternal(H,G,A){
  casadi_warning("OOQP interface still expermental");

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
  } else {
    // Split A in equalities and inequalities
    A_eq_.set(input(QP_A)(eq_,all_A_));
    A_ineq_.set(input(QP_A)(ineq_,all_A_));
    
    BA_eq_.set(input(QP_LBA)(eq_,all_A_));
    
    LBA_ineq_.set(input(QP_LBA)(ineq_,all_A_));
    UBA_ineq_.set(input(QP_UBA)(ineq_,all_A_));
    
    LBX_.set(input(QP_LBX));
    UBX_.set(input(QP_UBX));
    
    for (int k=0; k<LBA_ineq_.size();++k) 
      LBA_ineq_.data()[k] = iclow_[k]? LBA_ineq_.data()[k] : 0;
    for (int k=0; k<UBA_ineq_.size();++k) 
      UBA_ineq_.data()[k] = icupp_[k]? UBA_ineq_.data()[k] : 0;
    
    for (int k=0; k<LBX_.size();++k) 
      LBX_.data()[k] = ixlow_[k]? LBX_.data()[k] : 0;
    for (int k=0; k<UBX_.size();++k) 
      UBX_.data()[k] = ixupp_[k]? UBX_.data()[k] : 0;
    
    Hl_.set(input(QP_H)[Hl_nz_]);
  }
  
  // Just calling s_->solve repeatedly on a GondzioSolver is a memory leak
  // So we must always allocate a fresh solver
  if (s_) delete s_;
  s_ = new GondzioSolver( qp_, prob_ );
  
  s_->setMuTol(double(getOption("mutol")));
  s_->setArTol(double(getOption("artol")));
  
  int flag = s_->solve(prob_,vars_, resid_);
    
  vars_->x->copyIntoArray(&output(QP_PRIMAL).data()[0]);
  

  
  if (isnan(output(QP_PRIMAL).at(0))) {
    casadi_warning("nan in decision variables. You probably need to do a solver.reInit() call before evaluate().");
  }
  
  if(flag!=SUCCESSFUL_TERMINATION) ooqp_error("Solve",flag);

  
  output(QP_COST)[0] = prob_->objectiveValue(vars_);
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
  qp_ = new QpGenSparseMa27( nx, eq_.size(), ineq_.size(), input(QP_H).size() , input(QP_G).size(), input(QP_A).size() );
  

  // Split A in equalities and inequalities
  A_eq_   = input(QP_A)(eq_,all_A_);
  A_ineq_ = input(QP_A)(ineq_,all_A_);
  
  // Split up LBA & UBA in equalities and inequalities
  BA_eq_ = input(QP_LBA)(eq_,0);
  
  LBA_ineq_ = input(QP_LBA)(ineq_,0);
  UBA_ineq_ = input(QP_UBA)(ineq_,0);
  
  // Copy the decision variables bounds
  LBX_ = input(QP_LBX);
  UBX_ = input(QP_UBX);
  
  // Infinities on LBX & UBX are set to zero (required for OOQp)
  for (int k=0; k<input(QP_LBX).size(); ++k) {
    if (input(QP_LBX).at(k)==-numeric_limits<double>::infinity()) LBX_.at(k)=0;
  }
  for (int k=0; k<input(QP_UBX).size(); ++k) {
    if (input(QP_UBX).at(k)==numeric_limits<double>::infinity()) UBX_.at(k)=0;
  }
  
  
  // Set ixlow & ixupp: they have to be zero for infinite bounds and 1 otherwise
  for (int k=0; k<input(QP_LBX).size(); ++k) 
    ixlow_[k]=!(input(QP_LBX).at(k)==-numeric_limits<double>::infinity());
  for (int k=0; k<input(QP_UBX).size(); ++k) 
    ixupp_[k]= !(input(QP_UBX).at(k)==numeric_limits<double>::infinity());
  
  // Set iclow & icupp: they have to be zero for infinite bounds and 1 otherwise
  iclow_.resize(ineq_.size(),0);
  icupp_.resize(ineq_.size(),0);
  
  for (int k=0; k<LBA_ineq_.size();++k) 
    iclow_[k]=!(LBA_ineq_.at(k)==-numeric_limits<double>::infinity());
  
  for (int k=0; k<UBA_ineq_.size();++k) 
    icupp_[k]=!(UBA_ineq_.at(k)==numeric_limits<double>::infinity());

  // Infinities on LBA & UBA are set to zero (required for OOQp)
  for (int k=0; k<LBA_ineq_.size();k++) {
    if (LBA_ineq_.at(k)==-numeric_limits<double>::infinity()) LBA_ineq_.at(k)=0;
  } 
  
  for (int k=0; k<UBA_ineq_.size();k++) {
    if (UBA_ineq_.at(k)==numeric_limits<double>::infinity()) UBA_ineq_.at(k)=0;
  } 
  
  
  // Set Hl (the lower triangular part of H)
  Hl_.set(input(QP_H)[Hl_nz_]);
  
  // Because OOQP does not do const correctness properly, and because we cannot trust it, we copy all (rowind,col) data // NOTE: Joel: this appears unnecessary
  Hl_rowind_ = Hl_.sparsity().rowind();
  Hl_col_    = Hl_.sparsity().col();
  
  eq_rowind_ = A_eq_.sparsity().rowind();
  eq_col_ = A_eq_.sparsity().col();
  
  ineq_rowind_ = A_ineq_.sparsity().rowind();
  ineq_col_ = A_ineq_.sparsity().col();

  // Set all pointers to problem data
  prob_ = (QpGenData * )qp_->makeData( getPtr(input(QP_G)),
                                       getPtr(Hl_rowind_),  getPtr(Hl_col_),  getPtr(Hl_),
                                       getPtr(LBX_),  getPtr(ixlow_),
                                       getPtr(UBX_),  getPtr(ixupp_),
                                       getPtr(eq_rowind_), getPtr(eq_col_),  getPtr(A_eq_),
                                       getPtr(BA_eq_),
                                       getPtr(ineq_rowind_), getPtr(ineq_col_),  getPtr(A_ineq_),
                                       getPtr(LBA_ineq_),  getPtr(iclow_),
                                       getPtr(UBA_ineq_),  getPtr(icupp_));

  // Further setup of the QP problem
  vars_ = (QpGenVars *) qp_->makeVariables( prob_ );
  resid_ = (QpGenResiduals *) qp_->makeResiduals( prob_ );
 
  
}

void OOQPInternal::init(){
  
  QPSolverInternal::init();
  
  Hl_ = DMatrix(lowerSparsity(H));
  Hl_nz_ = lowerNZ(H);
  
  ixlow_.resize(nx,0);
  
  ixupp_.resize(nx,0);
  
  all_A_ = range(0,input(QP_A).size2());
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
  throw CasadiException(ss.str());
}

} // namespace Interfaces
} // namespace CasADi

