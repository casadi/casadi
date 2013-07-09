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

#include "sdp_sdqp_internal.hpp"

#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/fx/sx_function.hpp"

#include "symbolic/mx/mx_tools.hpp"
#include "symbolic/fx/mx_function.hpp"

using namespace std;
namespace CasADi {

SDPSDQPInternal* SDPSDQPInternal::clone() const{
  // Return a deep copy
  SDPSDQPInternal* node = new SDPSDQPInternal(st_);
  if(!node->is_init_)
    node->init();
  return node;
}
  
SDPSDQPInternal::SDPSDQPInternal(const std::vector<CRSSparsity> &st) : SDQPSolverInternal(st) {

  addOption("sdp_solver",       OT_SDPSOLVER, GenericType(), "The SDPSolver used to solve the SDQPs.");
  addOption("sdp_solver_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the SDPSOlver");
  
}

SDPSDQPInternal::~SDPSDQPInternal(){
}

void SDPSDQPInternal::evaluate(int nfdir, int nadir) {
  if (nfdir!=0 || nadir!=0) throw CasadiException("SDPSDQPInternal::evaluate() not implemented for forward or backward mode");
  
  cholesky_.setInput(input(SDQP_SOLVER_H));
  for (int k=0;k<cholesky_.input(0).size();++k) {
    cholesky_.input().at(k)*=0.5;
  }
  cholesky_.prepare();
  mapping_.setInput(cholesky_.getFactorization(true),"g_socp");
  std::copy(input(SDQP_SOLVER_C).begin(), input(SDQP_SOLVER_C).end(),mapping_.input("h_socp").begin());
  cholesky_.solveL(&mapping_.input("h_socp").data().front(),1,true);
  for (int k=0;k<mapping_.input("h_socp").size();++k) {
    mapping_.input("h_socp").at(k)*= 0.5;
  }
  
  mapping_.setInput(input(SDQP_SOLVER_F),"f_sdqp");
  mapping_.setInput(input(SDQP_SOLVER_G),"g_sdqp");
  
  mapping_.evaluate();
  
  std::copy(input(SDQP_SOLVER_A).begin(),input(SDQP_SOLVER_A).end(),sdpsolver_.input(SDP_SOLVER_A).begin());
  sdpsolver_.setInput(mapping_.output("f"),SDP_SOLVER_F);
  sdpsolver_.setInput(mapping_.output("g"),SDP_SOLVER_G);
  std::copy(input(SDQP_SOLVER_LBA).begin(),input(SDQP_SOLVER_LBA).end(),sdpsolver_.input(SDP_SOLVER_LBA).begin());
  std::copy(input(SDQP_SOLVER_UBA).begin(),input(SDQP_SOLVER_UBA).end(),sdpsolver_.input(SDP_SOLVER_UBA).begin());
  std::copy(input(SDQP_SOLVER_LBX).begin(),input(SDQP_SOLVER_LBX).end(),sdpsolver_.input(SDP_SOLVER_LBX).begin());
  std::copy(input(SDQP_SOLVER_UBX).begin(),input(SDQP_SOLVER_UBX).end(),sdpsolver_.input(SDP_SOLVER_UBX).begin());
  
  sdpsolver_.evaluate();
  
  std::copy(sdpsolver_.output(SDP_SOLVER_X).begin(),sdpsolver_.output(SDP_SOLVER_X).begin()+n_,output(SDQP_SOLVER_X).begin());
  setOutput(sdpsolver_.output(SDP_SOLVER_COST),SDQP_SOLVER_COST);
  setOutput(sdpsolver_.output(SDP_SOLVER_DUAL_COST),SDQP_SOLVER_DUAL_COST);
  if (!output(SDQP_SOLVER_DUAL).empty()) std::copy(sdpsolver_.output(SDP_SOLVER_DUAL).begin(),sdpsolver_.output(SDP_SOLVER_DUAL).begin()+output(SDQP_SOLVER_DUAL).size(),output(SDQP_SOLVER_DUAL).begin());
  if (!output(SDQP_SOLVER_P).empty()) std::copy(sdpsolver_.output(SDP_SOLVER_P).begin(),sdpsolver_.output(SDP_SOLVER_P).begin()+output(SDQP_SOLVER_P).size(),output(SDQP_SOLVER_P).begin());
  std::copy(sdpsolver_.output(SDP_SOLVER_X).begin(),sdpsolver_.output(SDP_SOLVER_X).begin()+n_,output(SDQP_SOLVER_X).begin());
  std::copy(sdpsolver_.output(SDP_SOLVER_LAM_A).begin(),sdpsolver_.output(SDP_SOLVER_LAM_A).end(),output(SDQP_SOLVER_LAM_A).begin());
  std::copy(sdpsolver_.output(SDP_SOLVER_LAM_X).begin(),sdpsolver_.output(SDP_SOLVER_LAM_X).begin()+n_,output(SDQP_SOLVER_LAM_X).begin());
}

void SDPSDQPInternal::init(){

  SDQPSolverInternal::init();
  cholesky_ = CSparseCholesky(st_[SDQP_STRUCT_H]);
  cholesky_.init();
  
  MX g_socp = msym("x", cholesky_.getFactorizationSparsity(true));
  MX h_socp = msym("h", n_);
  
  MX f_socp = sqrt(inner_prod(h_socp,h_socp));
  MX en_socp = 0.5/f_socp;

  MX f_sdqp = msym("f",input(SDQP_SOLVER_F).sparsity());
  MX g_sdqp = msym("g",input(SDQP_SOLVER_G).sparsity());
  
  std::vector<MX> fi(n_+1);
  MX znp = DMatrix(n_+1,n_+1);
  for (int k=0;k<n_;++k) {
    MX gk = vertcat(g_socp(ALL,k),DMatrix(1,1));
    MX fk = -blockcat(znp,gk,trans(gk),DMatrix(1,1));
    // TODO: replace with ALL
    fi.push_back(blkdiag(f_sdqp(range(f_sdqp.size2()*k,f_sdqp.size2()*(k+1)),range(f_sdqp.size2())),fk));
  }
  MX fin = en_socp*DMatrix::eye(n_+2);
  fin(n_,n_+1) = en_socp;
  fin(n_+1,n_) = en_socp;
  
  fi.push_back(blkdiag(DMatrix(f_sdqp.size2(),f_sdqp.size2()),-fin));
  
  MX h0 = vertcat(h_socp,DMatrix(1,1));
  MX g = blockcat(f_socp*DMatrix::eye(n_+1),h0,trans(h0),f_socp);
  
  g = blkdiag(g_sdqp,g);
  
  mapping_ = MXFunction(customIO("g_socp",g_socp,"h_socp",h_socp,"f_sdqp",f_sdqp,"g_sdqp",g_sdqp),customIO("f",vertcat(fi),"g",g));
  mapping_.init();

  // Create an sdpsolver instance
  SDPSolverCreator sdpsolver_creator = getOption("sdp_solver");
  sdpsolver_ = sdpsolver_creator(sdpStruct("g",mapping_.output("g").sparsity(),"f",mapping_.output("f").sparsity(),"a",horzcat(input(SDQP_SOLVER_A).sparsity(),sp_sparse(nc_,1))));

  if(hasSetOption("sdp_solver_options")){
    sdpsolver_.setOption(getOption("sdp_solver_options"));
  }
  
  // Initialize the SDP solver
  sdpsolver_.init();

  sdpsolver_.input(SDP_SOLVER_C).at(n_)=1;
  
  // Output arguments
  setNumOutputs(SDQP_SOLVER_NUM_OUT);
  output(SDQP_SOLVER_X) = DMatrix::zeros(n_,1);
  
  std::vector<int> r = range(input(SDQP_SOLVER_G).size1());
  output(SDQP_SOLVER_P) = sdpsolver_.output(SDP_SOLVER_P).empty() ? DMatrix() : sdpsolver_.output(SDP_SOLVER_P)(r,r);
  output(SDQP_SOLVER_DUAL) = sdpsolver_.output(SDP_SOLVER_DUAL).empty() ? DMatrix() : sdpsolver_.output(SDP_SOLVER_DUAL)(r,r);
  output(SDQP_SOLVER_COST) = 0.0;
  output(SDQP_SOLVER_DUAL_COST) = 0.0;
  output(SDQP_SOLVER_LAM_X) = DMatrix::zeros(n_,1);
  output(SDQP_SOLVER_LAM_A) = DMatrix::zeros(nc_,1);
 
}

} // namespace CasADi
