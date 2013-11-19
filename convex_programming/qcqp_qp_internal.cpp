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

#include "qcqp_qp_internal.hpp"

#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/fx/sx_function.hpp"

using namespace std;
namespace CasADi {

QCQPQPInternal* QCQPQPInternal::clone() const{
  // Return a deep copy
  QCQPQPInternal* node = new QCQPQPInternal(st_);
  if(!node->is_init_)
    node->init();
  return node;
}
  
QCQPQPInternal::QCQPQPInternal(const std::vector<CRSSparsity> &st) : QPSolverInternal(st) {

  addOption("qcqp_solver",       OT_QCQPSOLVER, GenericType(), "The QCQPSolver used to solve the QPs.");
  addOption("qcqp_solver_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the QCQPSOlver");
  
}

QCQPQPInternal::~QCQPQPInternal(){ 
}

void QCQPQPInternal::evaluate() {

  // Pass inputs of QP to QCQP form 
  qcqpsolver_.input(QCQP_SOLVER_A).set(input(QP_SOLVER_A));
  qcqpsolver_.input(QCQP_SOLVER_G).set(input(QP_SOLVER_G));
  qcqpsolver_.input(QCQP_SOLVER_H).set(input(QP_SOLVER_H));
  
  qcqpsolver_.input(QCQP_SOLVER_LBX).set(input(QP_SOLVER_LBX));
  qcqpsolver_.input(QCQP_SOLVER_UBX).set(input(QP_SOLVER_UBX));
  
  qcqpsolver_.input(QCQP_SOLVER_LBA).set(input(QP_SOLVER_LBA));
  qcqpsolver_.input(QCQP_SOLVER_UBA).set(input(QP_SOLVER_UBA));
  
  // Delegate computation to QCQP Solver
  qcqpsolver_.evaluate();
  
  // Pass the stats
  stats_["qcqp_solver_stats"] = qcqpsolver_.getStats();
  
  // Read the outputs from Ipopt
  output(QCQP_SOLVER_X).set(qcqpsolver_.output(QP_SOLVER_X));
  output(QCQP_SOLVER_COST).set(qcqpsolver_.output(QP_SOLVER_COST));
  output(QCQP_SOLVER_LAM_A).set(qcqpsolver_.output(QP_SOLVER_LAM_A));
  output(QCQP_SOLVER_LAM_X).set(qcqpsolver_.output(QP_SOLVER_LAM_X));
}

void QCQPQPInternal::init(){

  QPSolverInternal::init();

  // Create an qcqpsolver instance
  QCQPSolverCreator qcqpsolver_creator = getOption("qcqp_solver");
  qcqpsolver_ = qcqpsolver_creator(qcqpStruct("h",input(QP_SOLVER_H).sparsity(),"p",sp_sparse(0,n_),"a",input(QP_SOLVER_A).sparsity()));

  qcqpsolver_.setQPOptions();
  if(hasSetOption("qcqp_solver_options")){
    qcqpsolver_.setOption(getOption("qcqp_solver_options"));
  }
  
  // Initialize the NLP solver
  qcqpsolver_.init();
 
}

} // namespace CasADi

