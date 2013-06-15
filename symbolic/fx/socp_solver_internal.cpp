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

#include "socp_solver_internal.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../matrix/sparsity_tools.hpp"
#include "sx_function.hpp"
#include "../sx/sx_tools.hpp"
#include <numeric>

INPUTSCHEME(SOCPInput)
OUTPUTSCHEME(SOCPOutput)

using namespace std;
namespace CasADi{

// Constructor
SOCPSolverInternal::SOCPSolverInternal(const std::vector<CRSSparsity> &st) : st_(st) {
  addOption("ni",OT_INTEGERVECTOR, GenericType(), "Provide the size of each SOC constraint. Must sum up to N.");
  
  inputScheme_ = SCHEME_SOCPInput;
  outputScheme_ = SCHEME_SOCPOutput;

}
    
void SOCPSolverInternal::init() {
  // Call the init method of the base class
  FXInternal::init();
  
  ni_ = getOption("ni");
  m_ = ni_.size();
  
  const CRSSparsity& A = st_[SOCP_STRUCT_A];
  const CRSSparsity& G = st_[SOCP_STRUCT_G];
  
  N_ = std::accumulate(ni_.begin(), ni_.end(), 0);
  casadi_assert_message(N_==G.size1(),"SOCPSolverInternal: Supplied G sparsity: number of rows (" << G.size1() <<  ")  must match sum of vector provided with option 'ni' (" << N_ << ").");
  
  nc_ = A.size1();
  n_ = A.size2();
  
  casadi_assert_message(n_==G.size2(),"SOCPSolverInternal: Supplied G sparsity: number of cols (" << G.size2() <<  ")  must match number of decision variables (cols of A): " << n_ << ".");
  
  
  // Input arguments
  setNumInputs(SOCP_SOLVER_NUM_IN);
  input(SOCP_SOLVER_G) = DMatrix(G,0);
  input(SOCP_SOLVER_H) = DMatrix::zeros(N_,1);
  input(SOCP_SOLVER_E) = DMatrix::zeros(n_*m_,1);
  input(SOCP_SOLVER_F) = DMatrix::zeros(m_,1);
  input(SOCP_SOLVER_A) = DMatrix(A,0);
  input(SOCP_SOLVER_C) = DMatrix::zeros(n_);
  input(SOCP_SOLVER_LBX) = -DMatrix::inf(n_);
  input(SOCP_SOLVER_UBX) = DMatrix::inf(n_);
  input(SOCP_SOLVER_LBA) = -DMatrix::inf(nc_);
  input(SOCP_SOLVER_UBA) = DMatrix::inf(nc_);

  // Output arguments
  setNumOutputs(SOCP_SOLVER_NUM_OUT);
  output(SOCP_SOLVER_X) = DMatrix::zeros(n_,1);
  output(SOCP_SOLVER_COST) = 0.0;
  output(SOCP_SOLVER_LAM_X) = DMatrix::zeros(n_,1);
  output(SOCP_SOLVER_LAM_A) = DMatrix::zeros(nc_,1);
  
}

SOCPSolverInternal::~SOCPSolverInternal(){
}
 
void SOCPSolverInternal::evaluate(int nfdir, int nadir){
  throw CasadiException("SOCPSolverInternal::evaluate: Not implemented");
}
 
void SOCPSolverInternal::solve(){
  throw CasadiException("SOCPSolverInternal::solve: Not implemented");
}
 
} // namespace CasADi

  


