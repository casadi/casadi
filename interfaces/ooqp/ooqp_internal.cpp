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
#include "../../symbolic/fx/qp_solver.hpp"

#include "symbolic/matrix/sparsity_tools.hpp"
#include "symbolic/matrix/matrix_tools.hpp"
#include "symbolic/stl_vector_tools.hpp"

#include <Status.h>
// #define getPtr(x) (&(x).front())

using namespace std;
namespace CasADi {

/**
* Caveats:   lba = -inf && uba = inf is not allowed. Must be filtered out.
*
*/

  /**
   * \brief Reinitialize the problem 
   * This method needs to be called before evaluate() whenever the nature of any constraint has changed. This occurs when: \n
   *  - Any of LBA, UBA, LBX, UBX changes to/from (+-)infinity  \n
   *  - An entry of LBA becomes equal/unequal to UBA: this indicates that an inequality becomes an equality or visa versa. \n
   * 
   * You do not need to call this method before doing the very first evaluate() run
   */


/**
* The following transformation is performed on the problem:
*
*
* min f(x)
* s.t.  lbg <= g(x) <= ubg
*       lbx <=     x   <= ubx
*
* solve:
* min f([x,s])
*      s.t.  g(x) -s == 0
*      lbg <= s <= ubg
*      lbx <=     x   <= ubx
*/

OOQPInternal* OOQPInternal::clone() const{
  // Return a deep copy
  OOQPInternal* node = new OOQPInternal(st_);
  if(!node->is_init_)
    node->init();
  return node;
}
  
OOQPInternal::OOQPInternal(const std::vector<CRSSparsity>& st) : QPSolverInternal(st){
  addOption("print_level",OT_INTEGER,0,"Print level. OOQP listens to print_level 0, 10 and 100");
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

void OOQPInternal::evaluate() {
  if (inputs_check_) checkInputs();
  
  // Copy the bounds on X
  std::copy(input(QP_SOLVER_LBX).data().begin(),input(QP_SOLVER_LBX).data().end(),lbX_.begin());
  std::copy(input(QP_SOLVER_UBX).data().begin(),input(QP_SOLVER_UBX).data().end(),ubX_.begin());
  
  // Copy the constraint bounds
  std::copy(input(QP_SOLVER_LBA).data().begin(),input(QP_SOLVER_LBA).data().end(),lbX_.begin()+n_);
  std::copy(input(QP_SOLVER_UBA).data().begin(),input(QP_SOLVER_UBA).data().end(),ubX_.begin()+n_);
  
  // Infinities on LBX & UBX are set to zero (required for OOQp)
  // Set ixlow & ixupp: they have to be zero for infinite bounds and 1 otherwise
  for (int k=0; k<lbX_.size(); ++k) {
    ixlow_[k] = (lbX_[k] != -numeric_limits<double>::infinity());
    if (!ixlow_[k]) lbX_[k]=0;
    ixupp_[k] = (ubX_[k] != numeric_limits<double>::infinity());
    if (!ixupp_[k]) ubX_[k]=0;
  }
  
  // Obtain G
  std::copy(input(QP_SOLVER_G).data().begin(),input(QP_SOLVER_G).data().end(),G_.data().begin());
    
  // Obtain H
  input(QP_SOLVER_H).get(H_.data(),SPARSESYM);
  
  // Pass on QP_SOLVER_A
  vector<int> rowind,col;
  A_.sparsity().getSparsityCRS(rowind,col);
  int k_orig = 0;
  int k_new = 0;
  for(int r=0; r<rowind.size()-1; ++r) {
    for(int el=rowind[r]; el<rowind[r+1]; ++el){
      if (el<rowind[r+1]-1) {
        A_.data()[k_new] = input(QP_SOLVER_A).data()[k_orig];
        k_orig++;
      }
      k_new ++;
    }
  }

  // Get references to sparsity (NOTE: OOQP does not appear do be const correct)
  int *H_rowind = const_cast<int*>(getPtr(H_.sparsity().rowind()));
  int *H_col = const_cast<int*>(getPtr(H_.sparsity().col()));

  int *A_rowind = const_cast<int*>(getPtr(A_.sparsity().rowind()));
  int *A_col = const_cast<int*>(getPtr(A_.sparsity().col()));

  if (qp_) {
    delete qp_;
    delete prob_;
    delete vars_;
    delete resid_;
  }

  // Set up a Sparse solver; the decision space is [x,s]; there are no inequalities
  qp_ = new QpGenSparseMa27( n_ + nc_, nc_, 0, H_.size() , G_.size(), A_.size() );
  
  std::vector<int> rowind_empty_(1,0);
  
  // Set all pointers to problem data
  prob_ = (QpGenData * )qp_->makeData( getPtr(G_),
                                       H_rowind,  H_col,  getPtr(H_),
                                       getPtr(lbX_),  getPtr(ixlow_),
                                       getPtr(ubX_),  getPtr(ixupp_),
                                       A_rowind, A_col,  getPtr(A_),
                                       getPtr(b_),
                                       getPtr(rowind_empty_), 0,  0,
                                       0,  0,
                                       0,  0);
                                       
  // Further setup of the QP problem
  vars_ = (QpGenVars *) qp_->makeVariables( prob_ );

  resid_ = (QpGenResiduals *) qp_->makeResiduals( prob_ );
  
  // Temporary vector
  temp_.resize(n_+nc_,0);
  
  // Just calling s_->solve repeatedly on a GondzioSolver is a memory leak
  // So we must always allocate a fresh solver
  if (s_) delete s_;
  s_ = new GondzioSolver( qp_, prob_ );
  
  s_->setMuTol(double(getOption("mutol")));
  s_->setArTol(double(getOption("artol")));
  
  int flag = s_->solve(prob_,vars_, resid_);
    
  // Get primal solution
  vars_->x->copyIntoArray(getPtr(temp_));
  std::copy(temp_.begin(),temp_.begin()+n_,output(QP_SOLVER_X).begin());

  // Get multipliers for the bounds
//  vector<double> &lam_x = output(QP_SOLVER_LAM_X).data();
  
  // Set multipliers to zero
  output(QP_SOLVER_LAM_X).setAll(0);
  output(QP_SOLVER_LAM_A).setAll(0);
  
  // Lower bounds
  if (vars_->gamma->length() > 0) {
    vars_->gamma->copyIntoArray(getPtr(temp_));
    for (int k=0;k<n_;++k) 
      output(QP_SOLVER_LAM_X).data()[k] = -temp_[k];
    for (int k=0;k<nc_;++k) 
      output(QP_SOLVER_LAM_A).data()[k] = -temp_[n_+k];
  }

  // Upper bounds
  if (vars_->phi->length() > 0) {
    vars_->phi->copyIntoArray(getPtr(temp_));
    for (int k=0;k<n_;++k)
      output(QP_SOLVER_LAM_X).data()[k] += temp_[k];
    for (int k=0;k<nc_;++k) 
      output(QP_SOLVER_LAM_A).data()[k] += temp_[n_+k];
  }
  
  if(flag!=SUCCESSFUL_TERMINATION) ooqp_error("Solve",flag);
  
  output(QP_SOLVER_COST).at(0) = prob_->objectiveValue(vars_);

}

void OOQPInternal::init(){
   // Call the init method of the base class
  QPSolverInternal::init();
  
  // Reformulation of A
  std::vector<DMatrix> A_trans;
  A_trans.push_back(input(QP_SOLVER_A));
  A_trans.push_back(-DMatrix::eye(nc_));
  A_ = horzcat(A_trans);
  
  // Reformulation of G
  std::vector<DMatrix> G_trans;
  G_trans.push_back(input(QP_SOLVER_G));
  G_trans.push_back(DMatrix::zeros(nc_));
  G_ = vertcat(G_trans);
  
  // Reformulation of H
  std::vector<DMatrix> H_trans;
  H_trans.push_back(DMatrix(lowerSparsity(input(QP_SOLVER_H).sparsity())));
  H_trans.push_back(DMatrix(n_,nc_));
  
  H_ = horzcat(H_trans);
  H_trans.clear();
  H_trans.push_back(H_);
  H_trans.push_back(DMatrix(nc_,n_+nc_));
  
  H_ = vertcat(H_trans);
  
  // Allocate vectors to hold indicators of infinite variable bounds
  ixlow_.resize(n_+nc_,0);
  ixupp_.resize(n_+nc_,0);
  
  // Copy the decision variables bounds
  lbX_.resize(n_+nc_,0);
  ubX_.resize(n_+nc_,0);
  
  // The b vector of the transformation matrix is identcially zero
  b_.resize(n_+nc_,0);
  
  gOoqpPrintLevel = getOption("print_level");
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
  casadi_error(ss.str());
}

} // namespace CasADi

// #undef getPtr
