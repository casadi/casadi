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

#include "symbolic/mx/mx_tools.hpp"
#include "symbolic/fx/mx_function.hpp"

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

  // Create an MX for the decision variables
  MX X("X",nx_,1);
    
  // Put H, G, A sparsities in a vector...
  std::vector< CRSSparsity > sps;
  sps.push_back(input(QP_H).sparsity());
  sps.push_back(input(QP_G).sparsity());
  sps.push_back(input(QP_A).sparsity());
   
  // So that we can pass it on to createParent
  std::pair< MX, std::vector< MX > > mypair = createParent(sps);
  
  // V groups all parameters in an MX
  MX V(mypair.first);
  std::vector< MX > variables(mypair.second);
  
  // H_, G_, A_ depend on V
  H_=variables[0];
  G_=variables[1];
  A_=variables[2];

  // We're going to use two-argument objective and constraints to allow the use of parameters
  std::vector< MX > args;
  args.push_back(X);
  args.push_back(V);

  // The objective function looks exactly like a mathematical description of the NLP
  MXFunction QP_f(args, mul(trans(G_),X) + 0.5*mul(mul(trans(X),H_),X));
  QP_f.init();

  // So does the constraint function
  MXFunction QP_g(args, mul(A_,X));

  // Jacobian of the constraints
  MXFunction QP_j(args,A_);
  
/*  std:cout << (G_+mul(trans(H_),X)).dimString() << std::endl;*/
  // Gradient of the objective
  MXFunction QP_gf(args,G_+mul(H_,X));
  

  MX sigma("sigma");

  MX lambda("lambda",nc_,1);

  args.insert(args.begin()+1, lambda);
  args.insert(args.begin()+2, sigma);
  
  // Hessian of the Lagrangian
  MXFunction QP_h(args,H_*sigma);
  
  // Create an nlpsolver instance
  NLPSolverCreator nlpsolver_creator = getOption("nlp_solver");
  nlpsolver_ = nlpsolver_creator(QP_f,QP_g,QP_h,QP_j); // What to do with QP_gf?
  nlpsolver_.setQPOptions();
  if(hasSetOption("nlp_solver_options")){
    nlpsolver_.setOption(getOption("nlp_solver_options"));
  }
  nlpsolver_.setOption("parametric",true);
  
  // Initialize the NLP solver
  nlpsolver_.init();
 
}

} // namespace CasADi

