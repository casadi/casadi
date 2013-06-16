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

#include "socp_qcqp_internal.hpp"

#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/fx/sx_function.hpp"

using namespace std;
namespace CasADi {

SOCPQCQPInternal* SOCPQCQPInternal::clone() const{
  // Return a deep copy
  SOCPQCQPInternal* node = new SOCPQCQPInternal(st_);
  if(!node->is_init_)
    node->init();
  return node;
}
  
SOCPQCQPInternal::SOCPQCQPInternal(const std::vector<CRSSparsity> &st) : QCQPSolverInternal(st) {

  addOption("socp_solver",       OT_SOCPSOLVER, GenericType(), "The SOCPSolver used to solve the QCQPs.");
  addOption("socp_solver_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the SOCPSOlver");
  
}

SOCPQCQPInternal::~SOCPQCQPInternal(){
}

void SOCPQCQPInternal::evaluate(int nfdir, int nadir) {
  if (nfdir!=0 || nadir!=0) throw CasadiException("SOCPQCQPInternal::evaluate() not implemented for forward or backward mode");

  // Pass inputs of QCQP to SOCP form 
  socpsolver_.input(SOCP_SOLVER_A).set(input(QCQP_SOLVER_A));
  //socpsolver_.input(SOCP_SOLVER_G).set(input(QCQP_SOLVER_G));
  
  socpsolver_.input(SOCP_SOLVER_LBX).set(input(QCQP_SOLVER_LBX));
  socpsolver_.input(SOCP_SOLVER_UBX).set(input(QCQP_SOLVER_UBX));
  
  socpsolver_.input(SOCP_SOLVER_LBA).set(input(QCQP_SOLVER_LBA));
  socpsolver_.input(SOCP_SOLVER_UBA).set(input(QCQP_SOLVER_UBA));
  
  // Delegate computation to SOCP Solver
  socpsolver_.evaluate();
  
  // Read the outputs from Ipopt
  output(SOCP_SOLVER_X).set(socpsolver_.output(QCQP_SOLVER_X));
  output(SOCP_SOLVER_COST).set(socpsolver_.output(QCQP_SOLVER_COST));
  output(SOCP_SOLVER_LAM_A).set(socpsolver_.output(QCQP_SOLVER_LAM_A));
  output(SOCP_SOLVER_LAM_X).set(socpsolver_.output(QCQP_SOLVER_LAM_X));
}

void SOCPQCQPInternal::init(){

  QCQPSolverInternal::init();

  // Create an socpsolver instance
  SOCPSolverCreator socpsolver_creator = getOption("socp_solver");
  socpsolver_ = socpsolver_creator(socpStruct("g",sp_dense((nq_+1)*n_,n_),"a",input(QCQP_SOLVER_A).sparsity()));

  //socpsolver_.setQCQPOptions();
  if(hasSetOption("socp_solver_options")){
    socpsolver_.setOption(getOption("socp_solver_options"));
  }
  std::vector<int> ni(nq_+1);
  for (int i=0;i<nq_+1;++i) {
    ni[i] = n_;
  }
  socpsolver_.setOption("ni",ni);
  
  // Initialize the SOCP solver
  socpsolver_.init();
 
}

} // namespace CasADi
