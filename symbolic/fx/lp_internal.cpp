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

#include "lp_internal.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../matrix/sparsity_tools.hpp"

INPUTSCHEME(LPSolverInput)
OUTPUTSCHEME(LPSolverOutput)

using namespace std;
namespace CasADi{

// Constructor
LPSolverInternal::LPSolverInternal(const std::vector<CRSSparsity> &st) : st_(st) {
  casadi_assert_message(st_.size()==LP_STRUCT_NUM,"Problem structure mismatch");
  
  const CRSSparsity& A = st_[LP_STRUCT_A];
 
  n_ = A.size2();
  nc_ = A.size1();
  

  // Input arguments
  setNumInputs(LP_SOLVER_NUM_IN);
  input(LP_SOLVER_A) = DMatrix(A);
  input(LP_SOLVER_C) = DMatrix::zeros(n_);
  input(LP_SOLVER_LBA) = -DMatrix::inf(nc_);
  input(LP_SOLVER_UBA) = DMatrix::inf(nc_);
  input(LP_SOLVER_LBX) = -DMatrix::inf(n_);
  input(LP_SOLVER_UBX) = DMatrix::inf(n_);
  
  // Output arguments
  setNumOutputs(LP_SOLVER_NUM_OUT);
  output(LP_SOLVER_X) = DMatrix::zeros(n_);
  output(LP_SOLVER_COST) = 0.0;
  output(LP_SOLVER_LAM_X) = DMatrix::zeros(n_);
  output(LP_SOLVER_LAM_A) = DMatrix::zeros(nc_);
  
  inputScheme_ = SCHEME_LPSolverInput;
  outputScheme_ = SCHEME_LPSolverOutput;
}
    
void LPSolverInternal::init() {
  // Call the init method of the base class
  FXInternal::init();
}

LPSolverInternal::~LPSolverInternal(){
}
 
void LPSolverInternal::evaluate(int nfdir, int nadir){
  throw CasadiException("LPSolverInternal::evaluate: Not implemented");
}
 
void LPSolverInternal::solve(){
  throw CasadiException("LPSolverInternal::solve: Not implemented");
}
 
} // namespace CasADi

  


