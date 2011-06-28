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

#include "ooqp_internal.hpp"

#include "casadi/matrix/sparsity_tools.hpp"

using namespace std;
namespace CasADi {
namespace Interfaces {

OOQPInternal* OOQPInternal::clone() const{
  // Return a deep copy
  OOQPInternal* node = new OOQPInternal(H,G,A);
  if(!node->is_init)
    node->init();
  return node;
}
  
OOQPInternal::OOQPInternal(const CRSSparsity & H, const CRSSparsity & G, const CRSSparsity & A) : QPSolverInternal(H,G,A){
  addOption("constraints",OT_INTEGERVECTOR,std::vector<int>(nc,1),"Boolean vector indicating equality with 0 and inequality with 1");
  
  std::cout << "Warning: OOQP is highly experimental" << std::endl;
}

OOQPInternal::~OOQPInternal(){ 
}

void OOQPInternal::sort_constraints() {

}


void OOQPInternal::allocate() {
  // TODO: check if pointers are undefined
  
  
  
  if (hasOption("constraints")) {
    // Get constraints from Option
    constraints = getOption("constraints").toIntVector();
  } else {
    // Guess constraints
    for (int k=0; k<input(QP_LBX).size();k++) constraints[k] = !(input(QP_LBX).at(k)==input(QP_UBX).at(k));
  }
  
  
  // Find the number of inequalities
  n_ineq=0;
  for (int k=0;k<constraints.size();k++) n_ineq += constraints[k];
  
  // Find the number of equalities
  n_eq = constraints.size()-n_ineq;
  
  // Populate ineq
  ineq.resize(n_ineq);
  int cntineq=0;
  for (int k=0;k<constraints.size();k++) {
    if (constraints[k]) ineq[cntineq++]=k;
  }
  
  // Populate eq
  eq.resize(n_eq);
  int cnteq=0;
  for (int k=0;k<constraints.size();k++) {
    if (!constraints[k]) eq[cntineq++]=k;
  }
  
  qp = new QpGenSparseMa27( nx, n_eq, n_ineq, input(QP_H).size() , input(QP_G).size(), input(QP_A).size() );
  

  // Split A in equalities and inequalities
  A_eq   = input(QP_A)(eq,ALL);
  A_ineq = input(QP_A)(ineq,ALL);
  
  BA_eq  = input(QP_LBA)(eq,ALL);
  
  LBA_ineq = input(QP_LBA)(ineq,ALL);
  UBA_ineq = input(QP_UBA)(ineq,ALL);
  
  for (int k=0; k<input(QP_LBX).size();k++) ixlow[k]=!(input(QP_LBX).at(k)==-numeric_limits<double>::infinity());
  for (int k=0; k<input(QP_UBX).size();k++) ixupp[k]= !(input(QP_UBX).at(k)==numeric_limits<double>::infinity());
  
  iclow.resize(n_ineq,0);
  icupp.resize(n_ineq,0);
  
  for (int k=0; k<LBA_ineq.size();k++) iclow[k]=!(LBA_ineq.at(k)==-numeric_limits<double>::infinity());
  for (int k=0; k<UBA_ineq.size();k++) icupp[k]=!(UBA_ineq.at(k)==numeric_limits<double>::infinity());
  


  // TODO: makeData expects data pointers without a const modifier. So &Hl.sparsity().rowind()[0] is not allowed
  //       we have to make separate std::vector<int>
  
  //prob = (QpGenData * )qp->makeData( &input(QP_G).data()[0],
  //                       &Hl.sparsity().rowind()[0],  &Hl.sparsity().col()[0],  &Hl.data()[0],
  //                        &xlow[0],  &ixlow[0],
 //                        &xupp[0],  &ixupp[0],
  //                      &A_eq.sparsity().rowind()[0], &A_eq.sparsity().col()[0],  &A_eq.data()[0],
  //                     &BA_eq.data()[0],
  //                      &A_ineq.sparsity().rowind()[0], &A_ineq.sparsity().col()[0],  &A_ineq.data()[0],
  //                       &LBA_ineq.data()[0],  &iclow[0],
   //                      &UBA_ineq.data()[0],  &icupp[0]);

  vars = (QpGenVars *) qp->makeVariables( prob );
  resid = (QpGenResiduals *) qp->makeResiduals( prob );
  
  s     = new GondzioSolver( qp, prob );
  
  s->monitorSelf();
  int ierr = s->solve(prob,vars, resid);
  
}

void OOQPInternal::evaluate(int nfdir, int nadir) {

}

void OOQPInternal::init(){
  
  QPSolverInternal::init();
  
  if (hasSetOption("constraints")) {
    allocate();
  }
  
  Hl = DMatrix(lowerSparsity(H));
  
  xlow.resize(nx,0);
  ixlow.resize(nx,0);
  
  xupp.resize(nx,0);
  ixupp.resize(nx,0);
  
}

} // namespace Interfaces
} // namespace CasADi

