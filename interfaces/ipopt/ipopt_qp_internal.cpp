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

#include "ipopt_qp_internal.hpp"

#include "casadi/mx/mx_tools.hpp"
#include "casadi/fx/mx_function.hpp"

using namespace std;
namespace CasADi {
namespace Interfaces {

IpoptQPInternal* IpoptQPInternal::clone() const{
  // Return a deep copy
  IpoptQPInternal* node = new IpoptQPInternal(H,G,A);
  if(!node->is_init)
    node->init();
  return node;
}
  
IpoptQPInternal::IpoptQPInternal(const CRSSparsity & H, const CRSSparsity & G, const CRSSparsity & A) : QPSolverInternal(H,G,A){
  std::cout << "Warning: IPOPT QP is highly experimental" << std::endl;
}

IpoptQPInternal::~IpoptQPInternal(){ 
}

void IpoptQPInternal::evaluate(int nfdir, int nadir) {
  
  solver.input(NLP_P)[H_.mapping()] = input(QP_H);
  solver.input(NLP_P)[G_.mapping()] = input(QP_G);
  solver.input(NLP_P)[A_.mapping()] = input(QP_A);
  
  solver.input(NLP_LBX).set(input(QP_LBX));
  solver.input(NLP_UBX).set(input(QP_UBX));
  
  solver.input(NLP_LBG).set(input(QP_LBA));
  solver.input(NLP_UBG).set(input(QP_UBA));
  
  solver.evaluate();
  
  output(NLP_X_OPT).set(solver.output(NLP_X_OPT));
  output(NLP_COST).set(solver.output(NLP_COST));
}

void IpoptQPInternal::init(){
  
  QPSolverInternal::init();
  
  MX X("X",nx,1);
    
  std::vector< CRSSparsity > sps;
  sps.push_back(H);
  sps.push_back(G);
  sps.push_back(A);
  
  std::pair< MX, std::vector< MX > > mypair = createParent(sps);
  
  MX V(mypair.first);
  std::vector< MX > variables(mypair.second);
  
  H_=variables[0];
  G_=variables[1];
  A_=variables[2];
  
  
  std::vector< MX > args;
  args.push_back(X);
  args.push_back(V);
  
  MXFunction QP_f(args, prod(trans(G_),X) + 0.5*prod(prod(trans(X),H_),X));
  MXFunction QP_g(args, prod(A_,X));
  
  solver = IpoptSolver(QP_f,QP_g);

  solver.init();

}

} // namespace Interfaces
} // namespace CasADi

