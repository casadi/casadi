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

#include "sdp_socp_internal.hpp"

#include "casadi/symbolic/sx/sx_tools.hpp"
#include "casadi/symbolic/function/sx_function.hpp"
#include "casadi/symbolic/function/mx_function.hpp"
#include "casadi/symbolic/mx/mx_tools.hpp"

using namespace std;
namespace casadi {

SDPSOCPInternal* SDPSOCPInternal::clone() const{
  // Return a deep copy
  SDPSOCPInternal* node = new SDPSOCPInternal(st_);
  if(!node->is_init_)
    node->init();
  return node;
}

SDPSOCPInternal::SDPSOCPInternal(const std::vector<Sparsity> &st) : SOCPSolverInternal(st) {
  addOption("sdp_solver",       OT_SDPSOLVER, GenericType(),
            "The SDPSolver used to solve the SOCPs.");
  addOption("sdp_solver_options",       OT_DICTIONARY, GenericType(),
            "Options to be passed to the SDPSOlver");

}

SDPSOCPInternal::~SDPSOCPInternal(){
}

void SDPSOCPInternal::evaluate() {
  if (print_problem_) printProblem();

  mapping_.setInput(input(SOCP_SOLVER_G),0);
  mapping_.setInput(input(SOCP_SOLVER_H),1);
  mapping_.setInput(input(SOCP_SOLVER_E),2);
  mapping_.setInput(input(SOCP_SOLVER_F),3);

  mapping_.evaluate();

  sdpsolver_.setInput(mapping_.output(0),SDP_SOLVER_F);
  sdpsolver_.setInput(mapping_.output(1),SDP_SOLVER_G);
  sdpsolver_.setInput(input(SOCP_SOLVER_A),SDP_SOLVER_A);
  sdpsolver_.setInput(input(SOCP_SOLVER_C),SDP_SOLVER_C);
  sdpsolver_.setInput(input(SOCP_SOLVER_LBX),SDP_SOLVER_LBX);
  sdpsolver_.setInput(input(SOCP_SOLVER_UBX),SDP_SOLVER_UBX);
  sdpsolver_.setInput(input(SOCP_SOLVER_LBA),SDP_SOLVER_LBA);
  sdpsolver_.setInput(input(SOCP_SOLVER_UBA),SDP_SOLVER_UBA);

  sdpsolver_.evaluate();

  // Pass the stats
  stats_["sdp_solver_stats"] = sdpsolver_.getStats();

  setOutput(sdpsolver_.output(SDP_SOLVER_X),SOCP_SOLVER_X);
  setOutput(sdpsolver_.output(SDP_SOLVER_COST),SOCP_SOLVER_COST);
  setOutput(sdpsolver_.output(SDP_SOLVER_LAM_X),SOCP_SOLVER_LAM_X);
  setOutput(sdpsolver_.output(SDP_SOLVER_LAM_A),SOCP_SOLVER_LAM_A);
}

void SDPSOCPInternal::init(){

  SOCPSolverInternal::init();


  /*
  *   || Gi' x + hi ||_2 <=  ei'x + fi
  *
  *        <=>
  *
  *  | (ei' x + fi) I   Gi' x + hi   |   >= 0
  *  | (Gi' x + hi)'     ei' x + fi  |
  *
  *        <=>
  *
  *  | ei[k]I     Gi(:,k)'  |
  *  | Gi(:,k)     ei[k]  |
  *
  *      | k  (n)      \   i  (for each cone)
  *      v              \
  */


  MX G = MX::sym("G",input(SOCP_SOLVER_G).sparsity());
  MX H = MX::sym("H",input(SOCP_SOLVER_H).sparsity());
  MX E = MX::sym("E",input(SOCP_SOLVER_E).sparsity());
  MX F = MX::sym("F",input(SOCP_SOLVER_F).sparsity());


  int i_start;

  // The blocks that will make up Fi of SDP
  std::vector<MX> Fi;


  for (int k=0;k<n_;++k) {

    // The blocks that will end up on the diagonal of Fi of SDP
    std::vector<MX> Fi_d;

    i_start = 0;
    // Loop over all SOCP constraints
    for (int i=0;i<ni_.size();++i) {
      MX Gik = G(Slice(k),Slice(i_start,i_start+ni_[i])).T();
      MX Eik = E[n_*i+k];
      Fi_d.push_back(blockcat(Eik*MX::eye(ni_[i]),Gik,Gik.T(),Eik));
      i_start += ni_[i];
    }
    Fi.push_back(blkdiag(Fi_d));
  }
  // The blocks that will end up on the diagonal of G of SDP
  std::vector<MX> G_d;

  i_start = 0;
  // Loop over all SOCP constraints
  for (int i=0;i<ni_.size();++i) {
    MX Fi  = F[i];
    MX Hi  = H[range(i_start,i_start+ni_[i])];
    G_d.push_back(blockcat(Fi*MX::eye(ni_[i]),Hi,Hi.T(),Fi));
    i_start += ni_[i];
  }

  std::vector<MX> out;
  out.push_back(-horzcat(Fi));
  out.push_back(blkdiag(G_d));

  std::vector<MX> syms;
  syms.push_back(G);
  syms.push_back(H);
  syms.push_back(E);
  syms.push_back(F);

  mapping_ = MXFunction(syms,out);
  mapping_.init();

  log("SDPSOCPInternal::init","Created mapping function");

  // Create an sdpsolver instance
  SDPSolverCreator sdpsolver_creator = getOption("sdp_solver");
  sdpsolver_ = sdpsolver_creator(
    sdpStruct("a",input(SOCP_SOLVER_A).sparsity(),
              "f",mapping_.output(0).sparsity(),
              "g",mapping_.output(1).sparsity()));

  sdpsolver_.setSOCPOptions();
  if(hasSetOption("sdp_solver_options")){
    sdpsolver_.setOption(getOption("sdp_solver_options"));
  }

  // Initialize the SDP solver
  sdpsolver_.init();

  log("SDPSOCPInternal::init","Initialized SDP solver");
}

} // namespace casadi

