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

#include "sqic_internal.hpp"
#include "../../symbolic/fx/qp_solver.hpp"

#include "symbolic/matrix/sparsity_tools.hpp"
#include "symbolic/matrix/matrix_tools.hpp"
#include "symbolic/mx/mx_tools.hpp"
#include "symbolic/fx/mx_function.hpp"

#include "symbolic/stl_vector_tools.hpp"

#include "wsqic.hpp"

using namespace std;
namespace CasADi {

SQICInternal* SQICInternal::clone() const{
  // Return a deep copy
  SQICInternal* node = new SQICInternal(st_);
  if(!node->is_init_)
    node->init();
  return node;
}
  
SQICInternal::SQICInternal(const std::vector<CRSSparsity>& st) : QPSolverInternal(st){
  is_init_ = false;
}

SQICInternal::~SQICInternal(){ 
  sqicDestroy();
}

void SQICInternal::evaluate(int nfdir, int nadir) {
  casadi_assert_message(nfdir==0 && nadir==0, "SQIC::evaluate() not implemented for forward or backward mode");
  if (inputs_check_) checkInputs();

  std::copy(input(QP_SOLVER_X0).begin(),input(QP_SOLVER_X0).end(),x_.begin());
  
  std::transform(input(QP_SOLVER_LAM_X0).begin(),input(QP_SOLVER_LAM_X0).begin(),rc_.begin(),negate<double>());

  std::copy(input(QP_SOLVER_LBX).begin(),input(QP_SOLVER_LBX).end(),bl_.begin());
  std::copy(input(QP_SOLVER_UBX).begin(),input(QP_SOLVER_UBX).end(),bu_.begin());
 
  std::copy(input(QP_SOLVER_LBA).begin(),input(QP_SOLVER_LBA).end(),bl_.begin()+n_);
  std::copy(input(QP_SOLVER_UBA).begin(),input(QP_SOLVER_UBA).end(),bu_.begin()+n_);

  for (int i=0;i<n_+nc_+1;++i) {
    if (bl_[i]==-std::numeric_limits<double>::infinity()) bl_[i]=-inf_;
    if (bu_[i]==std::numeric_limits<double>::infinity()) bu_[i]=inf_;
  }
  
  formatA_.setInput(input(QP_SOLVER_A),0);
  formatA_.setInput(input(QP_SOLVER_G),1);
  formatA_.evaluate();
  
  sqicSolve(&output(QP_SOLVER_COST).data()[0]);

  std::copy(x_.begin(),x_.begin()+n_,output(QP_SOLVER_X).begin());
  std::transform(rc_.begin(),rc_.begin()+n_,output(QP_SOLVER_LAM_X).begin(),negate<double>());
  std::transform(rc_.begin()+n_,rc_.begin()+n_+nc_,output(QP_SOLVER_LAM_A).begin(),negate<double>());
  
  output(QP_SOLVER_COST)[0]+= x_[n_+nc_];
}

void SQICInternal::init(){
   // Call the init method of the base class
  QPSolverInternal::init();
  
  if (is_init_) sqicDestroy();
  
  inf_ = 1.0e+20;
  
  // Allocate data structures for SQIC
  bl_.resize(n_+nc_+1,0);
  bu_.resize(n_+nc_+1,0);
  x_.resize(n_+nc_+1,0);
  hs_.resize(n_+nc_+1,0);
  hEtype_.resize(n_+nc_+1,0);
  pi_.resize(nc_+1,0);
  rc_.resize(n_+nc_+1,0);
  
  locH_ = st_[QP_STRUCT_H].rowind();
  indH_ = st_[QP_STRUCT_H].col();
  
  // Fortran indices are one-based
  for (int i=0;i<indH_.size();++i) indH_[i]+=1;
  for (int i=0;i<locH_.size();++i) locH_[i]+=1;
  
  // Sparsity of augmented linear constraint matrix
  CRSSparsity A_ = trans(vertcat(st_[QP_STRUCT_A],sp_dense(1,n_)));
  locA_ = A_.rowind();
  indA_ = A_.col();
  
  // Fortran indices are one-based
  for (int i=0;i<indA_.size();++i) indA_[i]+=1;
  for (int i=0;i<locA_.size();++i) locA_[i]+=1;
  
  // helper functions for augmented linear constraint matrix
  MX a = msym("A",st_[QP_STRUCT_A]);
  MX g = msym("g",n_);
  std::vector<MX> ins;
  ins.push_back(a);
  ins.push_back(g);
  formatA_ = MXFunction(ins,trans(vertcat(a,trans(g))));
  formatA_.init();
  
  // Set objective row of augmented linear constraints
  bu_[n_+nc_] = inf_;
  bl_[n_+nc_] = -inf_;
  
  is_init_ = true;
  
  int n = n_;
  int m = nc_+1;
  
  int nnzA=formatA_.output().size();
  int nnzH=input(QP_SOLVER_H).size();
  
  std::fill(hEtype_.begin()+n_,hEtype_.end(),3);
    
  sqic(&m , &n, &nnzA, &indA_[0], &locA_[0], &formatA_.output().data()[0], &bl_[0], &bu_[0], &hEtype_[0], &hs_[0], &x_[0], &pi_[0], &rc_[0], &nnzH, &indH_[0], &locH_[0], &input(QP_SOLVER_H).data()[0]);
  
}

map<int,string> SQICInternal::calc_flagmap(){
  map<int,string> f;

  return f;
}
  
map<int,string> SQICInternal::flagmap = SQICInternal::calc_flagmap();

void SQICInternal::sqic_error(const string& module, int flag){
  // Find the error
  map<int,string>::const_iterator it = flagmap.find(flag);
  
  stringstream ss;
  if(it == flagmap.end()){
    ss << "Unknown error (" << flag << ") from module \"" << module << "\".";
  } else {
    ss << "Module \"" << module << "\" returned flag \"" << it->second << "\".";
  }
  ss << " Consult SQIC documentation.";
  casadi_error(ss.str());
}

} // namespace CasADi

// #undef getPtr
