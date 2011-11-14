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
#include "casadi/matrix/sparsity_tools.hpp"

INPUTSCHEME(QPInput)
OUTPUTSCHEME(QPOutput)

using namespace std;
namespace CasADi{

QPSolverInternal::QPSolverInternal() {
  //addOption("trans", OT_BOOLEAN, false);
}

// Constructor
QPSolverInternal::QPSolverInternal(const CRSSparsity &H, const CRSSparsity &A){
  addOption("convex", OT_BOOLEAN, false, "Specify true if you can guarantee that H will always be positive definite");

  nx_ = H.size2();
  nc_ = A.size1();
  
  casadi_assert_message(A.size2()==nx_ &&  H.size1()==H.size2(),
    "Got incompatible dimensions.   min          x'Hx + G'x s.t.   LBA <= Ax <= UBA :" << std::endl <<
    "H: " << H.dimString() << " - A: " << A.dimString() << std::endl <<
    "We need: H.size2()==A.size2(), H square & symmetric" << std::endl
  );

  // Sparsity
  CRSSparsity x_sparsity = sp_dense(nx_,1);
  CRSSparsity bounds_sparsity = sp_dense(nc_,1);
  
  // Input arguments
  setNumInputs(QP_NUM_IN);
  input(QP_X_INIT) = DMatrix(x_sparsity,0);
  input(QP_H) = DMatrix(H);
  input(QP_G) = DMatrix(x_sparsity);
  input(QP_A) = DMatrix(A);
  input(QP_LBA) = DMatrix(bounds_sparsity, -std::numeric_limits<double>::infinity());
  input(QP_UBA) = DMatrix(bounds_sparsity,  std::numeric_limits<double>::infinity());
  input(QP_LBX) = DMatrix(x_sparsity,      -std::numeric_limits<double>::infinity());
  input(QP_UBX) = DMatrix(x_sparsity,       std::numeric_limits<double>::infinity());
  
  // Output arguments
  setNumOutputs(QP_NUM_OUT);
  output(QP_PRIMAL) = DMatrix(x_sparsity);
  output(QP_COST) = 0.0;
  output(QP_DUAL_X) = DMatrix(x_sparsity);
  output(QP_DUAL_A) = DMatrix(bounds_sparsity);
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

  


