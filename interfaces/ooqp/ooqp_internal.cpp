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
#include "casadi/matrix/matrix_tools.hpp"

#include "casadi/pre_c99_support.hpp"

#include "Status.h"

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
  addOption("print_level",OT_INTEGER,0,"Print level. OOQP listends to print_level 0, 10 and 100");
  addOption("mutol",OT_REAL,1e-8,"tolerance as provided with setMuTol to OOQP");
  addOption("artol",OT_REAL,1e-8,"tolerance as provided with setArTol to OOQP");
  
  
  qp=0;
  prob=0;
  vars=0; 
  resid=0;
  s=0;
}

OOQPInternal::~OOQPInternal(){ 

  if (s) delete s;
  
  if (qp) {
    delete qp;
    delete prob;
    delete vars;
    delete resid;
  }
}

void OOQPInternal::evaluate(int nfdir, int nadir) {
  if (nfdir!=0 || nadir!=0) throw CasadiException("OOQPSolve::evaluate() not implemented for forward or backward mode");
  if (!qp) {
    allocate();
    assert(qp);
  } else {  
    // Split A in equalities and inequalities
    A_eq.set(input(QP_A)(eq,all_A));
    A_ineq.set(input(QP_A)(ineq,all_A));
    
    BA_eq.set(input(QP_LBA)(eq,all_A));
    
    LBA_ineq.set(input(QP_LBA)(ineq,all_A));
    UBA_ineq.set(input(QP_UBA)(ineq,all_A));
    
    LBX.set(input(QP_LBX));
    UBX.set(input(QP_UBX));
    
    for (int k=0; k<LBA_ineq.size();k++) LBA_ineq.data()[k] = iclow[k]? LBA_ineq.data()[k] : 0;
    for (int k=0; k<UBA_ineq.size();k++) UBA_ineq.data()[k] = icupp[k]? UBA_ineq.data()[k] : 0;
    
    for (int k=0; k<LBX.size();k++) LBX.data()[k] = ixlow[k]? LBX.data()[k] : 0;
    for (int k=0; k<UBX.size();k++) UBX.data()[k] = ixupp[k]? UBX.data()[k] : 0;
    
    Hl.set(input(QP_H)[Hl_nz]);
  }
  
  // Just calling s->solve repeatedly on a GondzioSolver is a memory leak
  // So we must always allocate a fresh solver
  if (s) delete s;
  s = new GondzioSolver( qp, prob );
  
  s->setMuTol(getOption("mutol").toDouble());
  s->setArTol(getOption("artol").toDouble());
  
  int flag = s->solve(prob,vars, resid);
    
  vars->x->copyIntoArray(&output(QP_X_OPT).data()[0]);
  

  
  if (isnan(output(QP_X_OPT).at(0))) {
     std::cerr << "WARNING: nan in decision variables. You probably need to do a solver.reInit() call before evaluate()." << std::endl;
  }
  
  if(flag!=SUCCESSFUL_TERMINATION) ooqp_error("Solve",flag);

  
  output(QP_COST)[0] = prob->objectiveValue(vars);
}

void OOQPInternal::allocate() {
  if (qp) {
    delete qp;
    delete prob;
    delete vars;
    delete resid;
  }
  
  
  // Decide if constraints are equality or inequality
  for (int k=0; k<input(QP_LBA).size();k++) constraints[k] = !(input(QP_LBA).at(k)==input(QP_UBA).at(k));
  
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
  
  
  // Set up a Sparse solver
  qp = new QpGenSparseMa27( nx, n_eq, n_ineq, input(QP_H).size() , input(QP_G).size(), input(QP_A).size() );
  

  // Split A in equalities and inequalities
  A_eq   = input(QP_A)(eq,all_A);
  A_ineq = input(QP_A)(ineq,all_A);
  
  // Split up LBA & UBA in equalities and inequalities
  BA_eq  = input(QP_LBA)(eq,0);
  
  LBA_ineq = input(QP_LBA)(ineq,0);
  UBA_ineq = input(QP_UBA)(ineq,0);
  
  // Copy the decision variables bounds
  LBX = input(QP_LBX);
  UBX = input(QP_UBX);

  
  // Infinities on LBX & UBX are set to zero (required for OOQp)
  for (int k=0; k<input(QP_LBX).size();k++) {
    if (input(QP_LBX).at(k)==-numeric_limits<double>::infinity()) LBX.at(k)=0;
  }
  for (int k=0; k<input(QP_UBX).size();k++) {
    if (input(QP_UBX).at(k)==numeric_limits<double>::infinity()) UBX.at(k)=0;
  }
  
  
  // Set ixlow & ixupp: they have to be zero for infinite bounds and 1 otherwise
  for (int k=0; k<input(QP_LBX).size();k++) ixlow[k]=!(input(QP_LBX).at(k)==-numeric_limits<double>::infinity());
  for (int k=0; k<input(QP_UBX).size();k++) ixupp[k]= !(input(QP_UBX).at(k)==numeric_limits<double>::infinity());
  

  // Set iclow & icupp: they have to be zero for infinite bounds and 1 otherwise
  iclow.resize(n_ineq,0);
  icupp.resize(n_ineq,0);
  
  for (int k=0; k<LBA_ineq.size();k++) iclow[k]=!(LBA_ineq.at(k)==-numeric_limits<double>::infinity());
  for (int k=0; k<UBA_ineq.size();k++) icupp[k]=!(UBA_ineq.at(k)==numeric_limits<double>::infinity());

  // Infinities on LBA & UBA are set to zero (required for OOQp)
  for (int k=0; k<LBA_ineq.size();k++) {
    if (LBA_ineq.at(k)==-numeric_limits<double>::infinity()) LBA_ineq.at(k)=0;
  } 
  
  for (int k=0; k<UBA_ineq.size();k++) {
    if (UBA_ineq.at(k)==numeric_limits<double>::infinity()) UBA_ineq.at(k)=0;
  } 
  
  
  // Set Hl (the lower triangular part of H)
  Hl.set(input(QP_H)[Hl_nz]);
  
  // Because OOQP does not do const correctness properly, and because we cannot trust it, we copy all (rowind,col) data
  Hl_rowind = Hl.sparsity().rowind();
  Hl_col    = Hl.sparsity().col();
  
  eq_rowind = A_eq.sparsity().rowind();
  eq_col = A_eq.sparsity().col();
  
  ineq_rowind = A_ineq.sparsity().rowind();
  ineq_col = A_ineq.sparsity().col();

  // Set all pointers to problem data
  prob = (QpGenData * )qp->makeData( &input(QP_G).data()[0],
                         vecptr(Hl_rowind),  vecptr(Hl_col),  &Hl.data()[0],
                          &LBX.data()[0],  vecptr(ixlow),
                         &UBX.data()[0],  vecptr(ixupp),
                       vecptr(eq_rowind), vecptr(eq_col),  &A_eq.data()[0],
                      &BA_eq.data()[0],
                       vecptr(ineq_rowind), vecptr(ineq_col),  &A_ineq.data()[0],
                        &LBA_ineq.data()[0],  vecptr(iclow),
                       &UBA_ineq.data()[0],  vecptr(icupp));

    
  // Further setup of the QP problem
  vars = (QpGenVars *) qp->makeVariables( prob );
  resid = (QpGenResiduals *) qp->makeResiduals( prob );
 
  
}

void OOQPInternal::init(){
  
  QPSolverInternal::init();
  
  Hl = DMatrix(lowerSparsity(H));
  Hl_nz = lowerNZ(H);
  
  ixlow.resize(nx,0);
  
  ixupp.resize(nx,0);
  
  all_A = range(0,input(QP_A).size2());
  
  constraints.resize(input(QP_A).size1(),0);
}

map<int,string> OOQPInternal::calc_flagmap(){
  map<int,string> f;

  f[SUCCESSFUL_TERMINATION] = "SUCCESSFUL_TERMINATION";
  f[NOT_FINISHED] = "NOT_FINISHED";
  f[MAX_ITS_EXCEEDED] = "MAX_ITS_EXCEEDED";
  f[INFEASIBLE] = "INFEASIBLE";
  f[UNKNOWN] = "UNKNOWN";
  return f;
}
  
map<int,string> OOQPInternal::flagmap = OOQPInternal::calc_flagmap();

void OOQPInternal::ooqp_error(const string& module, int flag){
  // Find the error
  map<int,string>::const_iterator it = flagmap.find(flag);
  
  stringstream ss;
  if(it == flagmap.end()){
    ss << "Unknown error (" << flag << ") from module \"" << module << "\".";
  } else {
    ss << "Module \"" << module << "\" returned flag \"" << it->second << "\".";
  }
  ss << " Consult OOQP documentation.";
  throw CasadiException(ss.str());
}

} // namespace Interfaces
} // namespace CasADi

