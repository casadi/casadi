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
#include "../matrix/matrix_tools.hpp"
#include "../matrix/sparsity_tools.hpp"

INPUTSCHEME(QPSolverInput)
OUTPUTSCHEME(QPSolverOutput)

using namespace std;
namespace CasADi{

// Constructor
QPSolverInternal::QPSolverInternal(const std::vector<CRSSparsity> &st) : st_(st) {

  casadi_assert_message(st_.size()==QP_STRUCT_NUM,"Problem structure mismatch");
  
  const CRSSparsity& A = st_[QP_STRUCT_A];
  const CRSSparsity& H = st_[QP_STRUCT_H];
  
  nx_ = H.size2();
  nc_ = A.isNull() ? 0 : A.size1();
  
  if (!A.isNull()) {
    casadi_assert_message(A.size2()==nx_,
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

  // Sparsity
  CRSSparsity x_sparsity = sp_dense(nx_,1);
  CRSSparsity bounds_sparsity = sp_dense(nc_,1);
  
  // Input arguments
  setNumInputs(QP_SOLVER_NUM_IN);
  input(QP_SOLVER_X0) = DMatrix(x_sparsity,0);
  input(QP_SOLVER_H) = DMatrix(H);
  input(QP_SOLVER_G) = DMatrix(x_sparsity);
  input(QP_SOLVER_A) = DMatrix(A);
  input(QP_SOLVER_LBA) = DMatrix(bounds_sparsity, -std::numeric_limits<double>::infinity());
  input(QP_SOLVER_UBA) = DMatrix(bounds_sparsity,  std::numeric_limits<double>::infinity());
  input(QP_SOLVER_LBX) = DMatrix(x_sparsity,      -std::numeric_limits<double>::infinity());
  input(QP_SOLVER_UBX) = DMatrix(x_sparsity,       std::numeric_limits<double>::infinity());
  
  // Output arguments
  setNumOutputs(QP_SOLVER_NUM_OUT);
  output(QP_SOLVER_X) = DMatrix(x_sparsity);
  output(QP_SOLVER_COST) = 0.0;
  output(QP_SOLVER_LAM_X) = DMatrix(x_sparsity);
  output(QP_SOLVER_LAM_A) = DMatrix(bounds_sparsity);
  
  inputScheme_ = SCHEME_QPSolverInput;
  outputScheme_ = SCHEME_QPSolverOutput;
}
    
void QPSolverInternal::init() {
  // Call the init method of the base class
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

  


