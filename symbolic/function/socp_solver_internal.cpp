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
SOCPSolverInternal::SOCPSolverInternal(const std::vector<Sparsity> &st) : st_(st) {
  addOption("ni",OT_INTEGERVECTOR, GenericType(), "Provide the size of each SOC constraint. Must sum up to N.");
  addOption("print_problem",OT_BOOLEAN,false,"Print out problem statement for debugging.");

  input_.scheme = SCHEME_SOCPInput;
  output_.scheme = SCHEME_SOCPOutput;

}
    
void SOCPSolverInternal::init() {
  // Call the init method of the base class
  FunctionInternal::init();
  
  ni_ = getOption("ni");
  print_problem_ = getOption("print_problem");
  
  m_ = ni_.size();
  
  const Sparsity& A = st_[SOCP_STRUCT_A];
  const Sparsity& G = st_[SOCP_STRUCT_G];
  
  N_ = std::accumulate(ni_.begin(), ni_.end(), 0);
  casadi_assert_message(N_==G.size2(),"SOCPSolverInternal: Supplied G sparsity: number of cols (" << G.size2() <<  ")  must match sum of vector provided with option 'ni' (" << N_ << ").");
  
  nc_ = A.size1();
  n_ = A.size2();
  
  casadi_assert_message(n_==G.size1(),"SOCPSolverInternal: Supplied G sparsity: number of rows (" << G.size1() <<  ")  must match number of decision variables (cols of A): " << n_ << ".");
  
  
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
 
void SOCPSolverInternal::evaluate(){
  throw CasadiException("SOCPSolverInternal::evaluate: Not implemented");
}
 
void SOCPSolverInternal::solve(){
  throw CasadiException("SOCPSolverInternal::solve: Not implemented");
}

void SOCPSolverInternal::printProblem(std::ostream &stream) const {
  stream << "SOCP Problem statement -- start" << std::endl;
  stream << "ni: "<< ni_ << std::endl;
  stream << "g: "<< std::endl;  input(SOCP_SOLVER_G).printDense(stream);
  stream << "h: "<< std::endl;  input(SOCP_SOLVER_H).printDense(stream);
  stream << "e: "<< std::endl;  input(SOCP_SOLVER_E).printDense(stream);
  stream << "f: "<< std::endl;  input(SOCP_SOLVER_F).printDense(stream);
  stream << "c: " << input(SOCP_SOLVER_C) << std::endl;
  stream << "a: " << input(SOCP_SOLVER_A) << std::endl;
  stream << "lba: " << input(SOCP_SOLVER_LBA) << std::endl;
  stream << "uba: " << input(SOCP_SOLVER_UBA) << std::endl;
  stream << "lbx: " << input(SOCP_SOLVER_LBX) << std::endl;
  stream << "ubx: " << input(SOCP_SOLVER_UBX) << std::endl;
  
  stream << "SOCP Problem statement -- end" << std::endl;
}
 

void SOCPSolverInternal::checkInputs() const {
  for (int i=0;i<input(SOCP_SOLVER_LBX).size();++i) {
    casadi_assert_message(input(SOCP_SOLVER_LBX).at(i)<=input(SOCP_SOLVER_UBX).at(i),"LBX[i] <= UBX[i] was violated for i=" << i << ". Got LBX[i]=" << input(SOCP_SOLVER_LBX).at(i) << " and UBX[i]=" << input(SOCP_SOLVER_UBX).at(i));
  }
  for (int i=0;i<input(SOCP_SOLVER_LBA).size();++i) {
    casadi_assert_message(input(SOCP_SOLVER_LBA).at(i)<=input(SOCP_SOLVER_UBA).at(i),"LBA[i] <= UBA[i] was violated for i=" << i << ". Got LBA[i]=" << input(SOCP_SOLVER_LBA).at(i) << " and UBA[i]=" << input(SOCP_SOLVER_UBA).at(i));
  }
}
 
} // namespace CasADi

  


