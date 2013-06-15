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

#include "qp_lp_internal.hpp"

#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/fx/sx_function.hpp"

using namespace std;
namespace CasADi {

QPLPInternal* QPLPInternal::clone() const{
  // Return a deep copy
  QPLPInternal* node = new QPLPInternal(st_);
  if(!node->is_init_)
    node->init();
  return node;
}
  
QPLPInternal::QPLPInternal(const std::vector<CRSSparsity> &st) : LPSolverInternal(st) {

  addOption("qp_solver",       OT_QPSOLVER, GenericType(), "The QPSOlver used to solve the LPs.");
  addOption("qp_solver_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the QPSOlver");
  
}

QPLPInternal::~QPLPInternal(){ 
}

void QPLPInternal::evaluate(int nfdir, int nadir) {
  if (nfdir!=0 || nadir!=0) throw CasadiException("QPLPInternal::evaluate() not implemented for forward or backward mode");

  // Pass inputs of LP to QP form 
  qpsolver_.input(QP_SOLVER_A).set(input(LP_SOLVER_A));
  qpsolver_.input(QP_SOLVER_G).set(input(LP_SOLVER_C));
  
  qpsolver_.input(QP_SOLVER_LBX).set(input(LP_SOLVER_LBX));
  qpsolver_.input(QP_SOLVER_UBX).set(input(LP_SOLVER_UBX));
  
  qpsolver_.input(QP_SOLVER_LBA).set(input(LP_SOLVER_LBA));
  qpsolver_.input(QP_SOLVER_UBA).set(input(LP_SOLVER_UBA));
  
  // Delegate computation to NLP Solver
  qpsolver_.evaluate();
  
  // Read the outputs from Ipopt
  output(QP_SOLVER_X).set(qpsolver_.output(LP_SOLVER_X));
  output(QP_SOLVER_COST).set(qpsolver_.output(LP_SOLVER_COST));
  output(QP_SOLVER_LAM_A).set(qpsolver_.output(LP_SOLVER_LAM_A));
  output(QP_SOLVER_LAM_X).set(qpsolver_.output(LP_SOLVER_LAM_X));
}

void QPLPInternal::init(){

  LPSolverInternal::init();

  // Create an qpsolver instance
  QPSolverCreator qpsolver_creator = getOption("qp_solver");
  qpsolver_ = qpsolver_creator(qpStruct("h",sp_sparse(n_,n_),"a",input(LP_SOLVER_A).sparsity()));

  qpsolver_.setLPOptions();
  if(hasSetOption("qp_solver_options")){
    qpsolver_.setOption(getOption("qp_solver_options"));
  }
  
  // Initialize the NLP solver
  qpsolver_.init();
 
}

} // namespace CasADi

