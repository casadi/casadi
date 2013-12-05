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

#include "stabilized_qp_solver_internal.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../matrix/sparsity_tools.hpp"

INPUTSCHEME(StabilizedQPSolverInput)
OUTPUTSCHEME(QPSolverOutput)

using namespace std;
namespace CasADi{

// Constructor
StabilizedQPSolverInternal::StabilizedQPSolverInternal(const std::vector<CRSSparsity> &st) : st_(st) {

  casadi_assert_message(st_.size()==QP_STRUCT_NUM,"Problem structure mismatch");
  
  const CRSSparsity& A = st_[QP_STRUCT_A];
  const CRSSparsity& H = st_[QP_STRUCT_H];
  
  n_ = H.size2();
  nc_ = A.isNull() ? 0 : A.size1();
  
  if (!A.isNull()) {
    casadi_assert_message(A.size2()==n_,
      "Got incompatible dimensions.   min          x'Hx + G'x s.t.   LBA <= Ax <= UBA :" << std::endl <<
      "H: " << H.dimString() << " - A: " << A.dimString() << std::endl <<
      "We need: H.size2()==A.size2()" << std::endl
    );
  } 
  
  casadi_assert_message(H==trans(H),
    "Got incompatible dimensions.   min          x'Hx + G'x" << std::endl <<
    "H: " << H.dimString() <<
    "We need H square & symmetric" << std::endl
  );

  // Sparsity
  CRSSparsity x_sparsity = sp_dense(n_,1);
  CRSSparsity bounds_sparsity = sp_dense(nc_,1);
  
  // Input arguments
  setNumInputs(STABILIZED_QP_SOLVER_NUM_IN);
  input(STABILIZED_QP_SOLVER_X0) = DMatrix(x_sparsity,0);
  input(STABILIZED_QP_SOLVER_H) = DMatrix(H);
  input(STABILIZED_QP_SOLVER_G) = DMatrix(x_sparsity);
  input(STABILIZED_QP_SOLVER_A) = DMatrix(A);
  input(STABILIZED_QP_SOLVER_LBA) = DMatrix(bounds_sparsity, -std::numeric_limits<double>::infinity());
  input(STABILIZED_QP_SOLVER_UBA) = DMatrix(bounds_sparsity,  std::numeric_limits<double>::infinity());
  input(STABILIZED_QP_SOLVER_LBX) = DMatrix(x_sparsity,      -std::numeric_limits<double>::infinity());
  input(STABILIZED_QP_SOLVER_UBX) = DMatrix(x_sparsity,       std::numeric_limits<double>::infinity());
  input(STABILIZED_QP_SOLVER_MUR) = DMatrix::zeros(1,1);
  input(STABILIZED_QP_SOLVER_MUE) = DMatrix(bounds_sparsity,0);
  input(STABILIZED_QP_SOLVER_MU) = DMatrix(bounds_sparsity,0);
  
  // Output arguments
  setNumOutputs(QP_SOLVER_NUM_OUT);
  output(QP_SOLVER_X) = DMatrix(x_sparsity);
  output(QP_SOLVER_COST) = 0.0;
  output(QP_SOLVER_LAM_X) = DMatrix(x_sparsity);
  output(QP_SOLVER_LAM_A) = DMatrix(bounds_sparsity);
  
  input_.scheme = SCHEME_StabilizedQPSolverInput;
  output_.scheme = SCHEME_QPSolverOutput;
}
    
void StabilizedQPSolverInternal::init() {
  // Call the init method of the base class
  FXInternal::init();
}

StabilizedQPSolverInternal::~StabilizedQPSolverInternal(){
}
 
void StabilizedQPSolverInternal::evaluate(){
  throw CasadiException("StabilizedQPSolverInternal::evaluate: Not implemented");
}
 
void StabilizedQPSolverInternal::solve(){
  throw CasadiException("StabilizedQPSolverInternal::solve: Not implemented");
}

void StabilizedQPSolverInternal::checkInputs() const {
  for (int i=0;i<input(STABILIZED_QP_SOLVER_LBX).size();++i) {
    casadi_assert_message(input(STABILIZED_QP_SOLVER_LBX).at(i)<=input(STABILIZED_QP_SOLVER_UBX).at(i),"LBX[i] <= UBX[i] was violated for i=" << i << ". Got LBX[i]=" << input(STABILIZED_QP_SOLVER_LBX).at(i) << " and UBX[i]=" << input(STABILIZED_QP_SOLVER_UBX).at(i));
  }
  for (int i=0;i<input(STABILIZED_QP_SOLVER_LBA).size();++i) {
    casadi_assert_message(input(STABILIZED_QP_SOLVER_LBA).at(i)<=input(STABILIZED_QP_SOLVER_UBA).at(i),"LBA[i] <= UBA[i] was violated for i=" << i << ". Got LBA[i]=" << input(STABILIZED_QP_SOLVER_LBA).at(i) << " and UBA[i]=" << input(STABILIZED_QP_SOLVER_UBA).at(i));
  }
}
 
} // namespace CasADi

  


