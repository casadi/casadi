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

#include "../../casadi/stl_vector_tools.hpp"
#include "../../casadi/matrix/matrix_tools.hpp"

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
  addOption("printLevel", OT_STRING,      GenericType(), "Defines the amount of text output during QP solution, see Section 5.7.: \"none\", \"low\", \"medium\" or \"high\"");
  
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
  qp_ = new qpOASES::SQProblem(nx_,nc_);
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
    flag = qp_->init(h,g,a,lb,ub,lbA,ubA,nWSR,cputime_ptr);
    called_once_ = true;
    casadi_assert(flag==qpOASES::SUCCESSFUL_RETURN || flag==qpOASES::RET_MAX_NWSR_REACHED);
  } else {
    flag = qp_->hotstart(h,g,a,lb,ub,lbA,ubA,nWSR, cputime_ptr);
    casadi_assert(flag==qpOASES::SUCCESSFUL_RETURN || flag==qpOASES::RET_MAX_NWSR_REACHED);
  }

  // Get optimal cost
  output(QP_COST).set(qp_->getObjVal());

  // Get the primal solution
  qp_->getPrimalSolution(&output(QP_PRIMAL).front());
  
  // Get the dual solution
  qp_->getDualSolution(&dual_.front());
  
  // Split up the dual solution in multipliers for the simple bounds and the linear bounds
  copy(dual_.begin(),   dual_.begin()+nx_,output(QP_DUAL_X).begin());
  copy(dual_.begin()+nx_,dual_.end(),     output(QP_DUAL_A).begin());
}

map<int,string> QPOasesInternal::calc_flagmap(){
  map<int,string> f;

  f[qpOASES::SUCCESSFUL_RETURN] = "SUCCESSFUL_RETURN";
  //f[NOT_FINISHED] = "NOT_FINISHED";
  f[qpOASES::RET_MAX_NWSR_REACHED] = "RET_MAX_NWSR_REACHED";
  f[qpOASES::RET_INIT_FAILED] = "RET INIT FAILED";
  f[qpOASES::RET_HOTSTART_FAILED] = "RET_HOTSTART_FAILED";
  return f;
}
  
map<int,string> QPOasesInternal::flagmap = QPOasesInternal::calc_flagmap();

void QPOasesInternal::qpoases_error(const string& module, int flag){
  // Find the error
  map<int,string>::const_iterator it = flagmap.find(flag);
  
  stringstream ss;
  if(it == flagmap.end()){
    ss << "Unknown error (" << flag << ") from module \"" << module << "\".";
  } else {
    ss << "Module \"" << module << "\" returned flag \"" << it->second << "\".";
  }
  ss << " Consult qpOASES documentation.";
  throw CasadiException(ss.str());
}

} // namespace Interfaces
} // namespace CasADi

