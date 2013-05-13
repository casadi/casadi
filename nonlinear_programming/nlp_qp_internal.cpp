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

#include "nlp_qp_internal.hpp"

#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/fx/sx_function.hpp"

using namespace std;
namespace CasADi {

NLPQPInternal* NLPQPInternal::clone() const{
  // Return a deep copy
  NLPQPInternal* node = new NLPQPInternal(input(QP_H).sparsity(),input(QP_A).sparsity());
  if(!node->is_init_)
    node->init();
  return node;
}
  
NLPQPInternal::NLPQPInternal(const CRSSparsity& H_, const CRSSparsity &A_) : QPSolverInternal(H_,A_) {

  addOption("nlp_solver",       OT_NLPSOLVER, GenericType(), "The NLPSOlver used to solve the QPs.");
  addOption("nlp_solver_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the NLPSOlver");
  
}

NLPQPInternal::~NLPQPInternal(){ 
}

void NLPQPInternal::evaluate(int nfdir, int nadir) {
  if (nfdir!=0 || nadir!=0) throw CasadiException("NLPQPInternal::evaluate() not implemented for forward or backward mode");

  int k = 0;
  
 // Pass inputs of QP to NLP form 
  
  std::copy(input(QP_H).data().begin(),input(QP_H).data().end(),nlpsolver_.input(NLP_SOLVER_P).data().begin()+k); k+= input(QP_H).size();
  std::copy(input(QP_G).data().begin(),input(QP_G).data().end(),nlpsolver_.input(NLP_SOLVER_P).data().begin()+k); k+= input(QP_G).size();
  std::copy(input(QP_A).data().begin(),input(QP_A).data().end(),nlpsolver_.input(NLP_SOLVER_P).data().begin()+k);
  

  nlpsolver_.input(NLP_SOLVER_LBX).set(input(QP_LBX));
  nlpsolver_.input(NLP_SOLVER_UBX).set(input(QP_UBX));
  
  nlpsolver_.input(NLP_SOLVER_LBG).set(input(QP_LBA));
  nlpsolver_.input(NLP_SOLVER_UBG).set(input(QP_UBA));
  
  // Delegate computation to NLP Solver
  nlpsolver_.evaluate();
  
  // Read the outputs from Ipopt
  output(QP_PRIMAL).set(nlpsolver_.output(NLP_SOLVER_X));
  output(QP_COST).set(nlpsolver_.output(NLP_SOLVER_F));
  output(QP_LAMBDA_A).set(nlpsolver_.output(NLP_SOLVER_LAM_G));
  output(QP_LAMBDA_X).set(nlpsolver_.output(NLP_SOLVER_LAM_X));
}

void NLPQPInternal::init(){

  
  QPSolverInternal::init();

  // Create a symbolic matrix for the decision variables
  SXMatrix X = ssym("X",nx_,1);

  // Parameters to the problem
  SXMatrix H = ssym("H",input(QP_H).sparsity());
  SXMatrix G = ssym("G",input(QP_G).sparsity());
  SXMatrix A = ssym("A",input(QP_A).sparsity());

  // Put parameters in a vector
  std::vector< SXMatrix > par;
  par.push_back(H.data());
  par.push_back(G.data());
  par.push_back(A.data());

  // The nlp looks exactly like a mathematical description of the NLP
  SXFunction QP_nlp(nlpIn("x",X,"p",vertcat(par)),
		    nlpOut("f",mul(trans(G),X) + 0.5*mul(mul(trans(X),H),X),
			   "g",mul(A,X)));

  // Create an nlpsolver instance
  NLPSolverCreator nlpsolver_creator = getOption("nlp_solver");
  nlpsolver_ = nlpsolver_creator(QP_nlp);

  nlpsolver_.setQPOptions();
  if(hasSetOption("nlp_solver_options")){
    nlpsolver_.setOption(getOption("nlp_solver_options"));
  }
  
  // Initialize the NLP solver
  nlpsolver_.init();
 
}

} // namespace CasADi

