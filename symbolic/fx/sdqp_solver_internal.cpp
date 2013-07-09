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

#include "sdqp_solver_internal.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../matrix/sparsity_tools.hpp"
#include "sx_function.hpp"
#include "../sx/sx_tools.hpp"

INPUTSCHEME(SDQPInput)
OUTPUTSCHEME(SDQPOutput)

using namespace std;
namespace CasADi{

// Constructor
SDQPSolverInternal::SDQPSolverInternal(const std::vector<CRSSparsity> &st) : st_(st) {

  addOption("sdp_solver",       OT_SDPSOLVER, GenericType(), "The SDQPSolver used to solve the SDPs.");
  addOption("sdp_solver_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the SDPSOlver");
  
  casadi_assert_message(st_.size()==SDQP_STRUCT_NUM,"Problem structure mismatch");
  
  const CRSSparsity& A = st_[SDQP_STRUCT_A];
  const CRSSparsity& G = st_[SDQP_STRUCT_G];
  const CRSSparsity& F = st_[SDQP_STRUCT_F];
  const CRSSparsity& H = st_[SDQP_STRUCT_H];
  
  casadi_assert_message(G==G.transpose(),"SDQPSolverInternal: Supplied G sparsity must symmetric but got " << G.dimString());
  casadi_assert_message(H==H.transpose(),"SDQPSolverInternal: Supplied H sparsity must symmetric but got " << H.dimString());
  
  m_ = G.size1();
  
  nc_ = A.size1();
  n_ = H.size1();
  
  casadi_assert_message(F.size2()==m_,"SDQPSolverInternal: Supplied F sparsity: number of columns (" << F.size2() <<  ")  must match m (" << m_ << ")");
  
  casadi_assert_message(A.size2()==n_,"SDQPSolverInternal: Supplied A sparsity: number of columns (" << A.size2() <<  ")  must match n (" << n_ << ")");
  
  casadi_assert_message(F.size1()%n_==0,"SDQPSolverInternal: Supplied F sparsity: number of rows (" << F.size2() <<  ")  must be an integer multiple of n (" << n_ << "), but got remainder " << F.size1()%n_);
  
  // Input arguments
  setNumInputs(SDQP_SOLVER_NUM_IN);
  input(SDQP_SOLVER_H) = DMatrix(H,0);
  input(SDQP_SOLVER_G) = DMatrix(G,0);
  input(SDQP_SOLVER_F) = DMatrix(F,0);
  input(SDQP_SOLVER_A) = DMatrix(A,0);
  input(SDQP_SOLVER_C) = DMatrix::zeros(n_);
  input(SDQP_SOLVER_LBX) = -DMatrix::inf(n_);
  input(SDQP_SOLVER_UBX) = DMatrix::inf(n_);
  input(SDQP_SOLVER_LBA) = -DMatrix::inf(nc_);
  input(SDQP_SOLVER_UBA) = DMatrix::inf(nc_);

  for (int i=0;i<n_;i++) {
    CRSSparsity s = input(SDQP_SOLVER_F)(range(i*m_,(i+1)*m_),ALL).sparsity();
    casadi_assert_message(s==s.transpose(),"SDQPSolverInternal: Each supplied Fi must be symmetric. But got " << s.dimString() <<  " for i = " << i << ".");
  }
  
  inputScheme_ = SCHEME_SDQPInput;
  outputScheme_ = SCHEME_SDQPOutput;

}
    
void SDQPSolverInternal::init() {
  // Call the init method of the base class
  FXInternal::init();
  
}

SDQPSolverInternal::~SDQPSolverInternal(){
}

void SDQPSolverInternal::printProblem(std::ostream &stream) const {
  stream << "SDQP Problem statement -- start" << std::endl;
  
  stream << "h: "<< std::endl;  input(SDQP_SOLVER_H).printDense(stream);
  stream << "c: "<< std::endl;  input(SDQP_SOLVER_C).printDense(stream);
  stream << "f: "<< std::endl;  input(SDQP_SOLVER_F).printDense(stream);
  stream << "g: "<< std::endl;  input(SDQP_SOLVER_G).printDense(stream);
  stream << "a: "<< std::endl;  input(SDQP_SOLVER_A).printDense(stream);
  stream << "lba: " << input(SDQP_SOLVER_LBA) << std::endl;
  stream << "uba: " << input(SDQP_SOLVER_UBA) << std::endl;
  stream << "lbx: " << input(SDQP_SOLVER_LBX) << std::endl;
  stream << "ubx: " << input(SDQP_SOLVER_UBX) << std::endl;
  
  stream << "SDQP Problem statement -- end" << std::endl;
}
 
void SDQPSolverInternal::evaluate(int nfdir, int nadir){
  throw CasadiException("SDQPSolverInternal::evaluate: Not implemented");
}
 
void SDQPSolverInternal::solve(){
  throw CasadiException("SDQPSolverInternal::solve: Not implemented");
}

} // namespace CasADi

  


