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

#include "qp_solver_internal.hpp"
#include "casadi/matrix/matrix_tools.hpp"

using namespace std;
namespace CasADi{

QPSolverInternal::QPSolverInternal() {
  //addOption("trans", OT_BOOLEAN, false);
}

// Constructor
QPSolverInternal::QPSolverInternal(const CRSSparsity &H_, const CRSSparsity &G_, const CRSSparsity &A_) : H(H_), G(G_), A(A_) {
  int nx = H.size2();
  if (G.size1()!=nx || A.size2()!=0 || H.size1()!=H.size2()) {
    stringstream ss;
    ss << "Got incompatible dimensions.   min          x'Hx + G'x s.t.   LBA <= Ax <= UBA :" << std::endl;
    ss << "H: " << H.dimString() << " - G: " << G.dimString() << " - A: " << A.dimString() << std::endl;
    throw CasadiException(ss.str());
  }
  if (G.numel() != G.size()) {
    stringstream ss;
    ss << "G in  min          x'Hx + G'x s.t.   LBA <= Ax <= UBA :" << std::endl;
    ss << "must be dense. Got " << G.dimString() <<  std::endl;
    throw CasadiException(ss.str());
  }

}
    
void QPSolverInternal::init() {
  setNumInputs(QP_NUM_IN);
  input_.resize(QP_NUM_IN);
  
  input(QP_X_INIT) = DMatrix(nx,1,0);
    
  input(QP_H) = DMatrix(H);
  input(QP_G) = DMatrix(G);
  input(QP_A) = DMatrix(A);
  
  input(QP_LBA) = DMatrix(dot(input(QP_X_INIT),input(QP_A)).sparsity());
  input(QP_UBA) = DMatrix(input(QP_LBA).sparsity());
  
  input(QP_LBX) = DMatrix(input(QP_X_INIT).sparsity());
  input(QP_UBX) = DMatrix(input(QP_X_INIT).sparsity());
  
  FXInternal::init();
}

QPSolverInternal::~QPSolverInternal(){
}
 
void QPSolverInternal::evaluate(int nfdir, int nadir){
  throw CasadiException("QPSolverInternal::evaluate: Not implemented");
}
 
void QPSolverInternal::solve(){
  throw CasadiException("QPSolverInternal::solve: Not implemented");
}
 
} // namespace CasADi

  


