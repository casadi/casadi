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

#include "qcqp_solver_internal.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../matrix/sparsity_tools.hpp"

INPUTSCHEME(QCQPSolverInput)
OUTPUTSCHEME(QCQPSolverOutput)

using namespace std;
namespace CasADi{

// Constructor
QCQPSolverInternal::QCQPSolverInternal(const std::vector<CRSSparsity> &st) : st_(st) {

  casadi_assert_message(st_.size()==QCQP_STRUCT_NUM,"Problem structure mismatch");
  
  const CRSSparsity& A = st_[QCQP_STRUCT_A];
  const CRSSparsity& P = st_[QCQP_STRUCT_P];
  const CRSSparsity& H = st_[QCQP_STRUCT_H];
  
  n_ = H.size2();
  nc_ = A.isNull() ? 0 : A.size1();
  
  if (!A.isNull()) {
    casadi_assert_message(A.size2()==n_,
      "Got incompatible dimensions.   min          x'Hx + G'x s.t.   LBA <= Ax <= UBA :" << std::endl <<
      "H: " << H.dimString() << " - A: " << A.dimString() << std::endl <<
      "We need: H.size2()==A.size2()" << std::endl
    );
  } 
  
  casadi_assert_message(H.size1()==H.size2(),
    "Got incompatible dimensions.   min          x'Hx + G'x" << std::endl <<
    "H: " << H.dimString() <<
    "We need H square & symmetric" << std::endl
  );
  
  casadi_assert_message(P.size2()==n_,"Got incompatible dimensions. Number of colums in P (" << P.size2() << ") must match n (" << n_ << ").");

  casadi_assert_message(P.size1() % n_ == 0,"Got incompatible dimensions. Number of rows in P (" << P.size1() << ") must be a multiple of n (" << n_ << ").");
  
  
  nq_ = P.size1() / n_;

  // Sparsity
  CRSSparsity x_sparsity = sp_dense(n_,1);
  CRSSparsity bounds_sparsity = sp_dense(nc_,1);
  
  // Input arguments
  setNumInputs(QCQP_SOLVER_NUM_IN);
  input(QCQP_SOLVER_X0) = DMatrix(x_sparsity,0);
  input(QCQP_SOLVER_H) = DMatrix(H);
  input(QCQP_SOLVER_G) = DMatrix(x_sparsity);
  input(QCQP_SOLVER_A) = DMatrix(A);
  input(QCQP_SOLVER_P) = DMatrix(P);
  input(QCQP_SOLVER_Q) = DMatrix::zeros(nq_*n_,1);
  input(QCQP_SOLVER_R) = DMatrix::zeros(nq_,1);  
  input(QCQP_SOLVER_LBA) = DMatrix(bounds_sparsity, -std::numeric_limits<double>::infinity());
  input(QCQP_SOLVER_UBA) = DMatrix(bounds_sparsity,  std::numeric_limits<double>::infinity());
  input(QCQP_SOLVER_LBX) = DMatrix(x_sparsity,      -std::numeric_limits<double>::infinity());
  input(QCQP_SOLVER_UBX) = DMatrix(x_sparsity,       std::numeric_limits<double>::infinity());

  for (int i=0;i<nq_;++i) {
    DMatrix Pi = input(QCQP_SOLVER_P)(range(i*n_,(i+1)*n_),ALL);
    casadi_assert_message(Pi.sparsity()==trans(Pi.sparsity()),"We need Pi square & symmetric but got " << Pi.dimString() << " for i = " << i << ".");
  }
  
  // Output arguments
  setNumOutputs(QCQP_SOLVER_NUM_OUT);
  output(QCQP_SOLVER_X) = DMatrix(x_sparsity);
  output(QCQP_SOLVER_COST) = 0.0;
  output(QCQP_SOLVER_LAM_X) = DMatrix(x_sparsity);
  output(QCQP_SOLVER_LAM_A) = DMatrix(bounds_sparsity);
  
  input_.scheme = SCHEME_QCQPSolverInput;
  output_.scheme = SCHEME_QCQPSolverOutput;
}
    
void QCQPSolverInternal::init() {
  // Call the init method of the base class
  FXInternal::init();
}

QCQPSolverInternal::~QCQPSolverInternal(){
}
 
void QCQPSolverInternal::evaluate(){
  throw CasadiException("QCQPSolverInternal::evaluate: Not implemented");
}
 
void QCQPSolverInternal::solve(){
  throw CasadiException("QCQPSolverInternal::solve: Not implemented");
}

void QCQPSolverInternal::checkInputs() const {
  for (int i=0;i<input(QCQP_SOLVER_LBX).size();++i) {
    casadi_assert_message(input(QCQP_SOLVER_LBX).at(i)<=input(QCQP_SOLVER_UBX).at(i),"LBX[i] <= UBX[i] was violated for i=" << i << ". Got LBX[i]=" << input(QCQP_SOLVER_LBX).at(i) << " and UBX[i]=" << input(QCQP_SOLVER_UBX).at(i));
  }
  for (int i=0;i<input(QCQP_SOLVER_LBA).size();++i) {
    casadi_assert_message(input(QCQP_SOLVER_LBA).at(i)<=input(QCQP_SOLVER_UBA).at(i),"LBA[i] <= UBA[i] was violated for i=" << i << ". Got LBA[i]=" << input(QCQP_SOLVER_LBA).at(i) << " and UBA[i]=" << input(QCQP_SOLVER_UBA).at(i));
  }
}
 
} // namespace CasADi

  


