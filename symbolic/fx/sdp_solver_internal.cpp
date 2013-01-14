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

#include "sdp_solver_internal.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../matrix/sparsity_tools.hpp"
#include "sx_function.hpp"
#include "../sx/sx_tools.hpp"

INPUTSCHEME(SDPInput)
OUTPUTSCHEME(SDPOutput)

using namespace std;
namespace CasADi{

SDPSolverInternal::SDPSolverInternal() {
  
  inputScheme = SCHEME_SDPInput;
  outputScheme = SCHEME_SDPOutput;
}



// Constructor
SDPSolverInternal::SDPSolverInternal(const CRSSparsity &C, const CRSSparsity &A) {
  addOption("calc_p",OT_BOOLEAN, true, "Indicate if the P-part of primal solution should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (n x n).");
  addOption("calc_dual",OT_BOOLEAN, true, "Indicate if dual should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (n x n).");
  
  casadi_assert_message(C==C.transpose(),"SDPSolverInternal: Supplied C sparsity must symmetric but got " << A.dimString());
  
  n_ = C.size1();
  
  casadi_assert_message(A.size2()==n_,"SDPSolverInternal: Supplied A sparsity: number of columns (" << A.size2() <<  ")  must match n (" << n_ << ")");
  
  casadi_assert_message(A.size1()%n_==0,"SDPSolverInternal: Supplied A sparsity: number of rows (" << A.size2() <<  ")  must be an integer multiple of n (" << n_ << "), but got remainder " << A.size1()%n_);
  
  m_ = A.size1()/n_;
  
  // Input arguments
  setNumInputs(SDP_NUM_IN);
  input(SDP_C) = DMatrix(C,0);
  input(SDP_A) = DMatrix(A,0);
  input(SDP_B) = DMatrix::zeros(m_);

  for (int i=0;i<m_;i++) {
    CRSSparsity s = input(SDP_A)(range(i*n_,(i+1)*n_),ALL).sparsity();
    casadi_assert_message(s==s.transpose(),"SDPSolverInternal: Each supplied Ai must be symmetric. But got " << s.dimString() <<  " for i = " << i << ".");
  }

}
    
void SDPSolverInternal::init() {
  // Call the init method of the base class
  FXInternal::init();
  
  calc_p_ = getOption("calc_p");
  calc_dual_ = getOption("calc_dual");
  
  // Output arguments
  setNumOutputs(SDP_NUM_OUT);
  output(SDP_PRIMAL) = DMatrix::zeros(m_,1);
  output(SDP_PRIMAL_P) = calc_p_? DMatrix::zeros(n_,n_) : DMatrix();
  output(SDP_DUAL) = calc_dual_? DMatrix::zeros(n_,n_) : DMatrix();
  output(SDP_PRIMAL_COST) = 0.0;
  output(SDP_DUAL_COST) = 0.0;
  
  SXMatrix var = ssym("var",input(SDP_A).sparsity());
  std::vector<SXMatrix> out(m_);
  for (int i=0;i<m_;++i) {
    out[i] = var(range(i*n_,(i+1)*n_),ALL);
  }
  mapping_ = SXFunction(var,out);
  mapping_.init();
}

SDPSolverInternal::~SDPSolverInternal(){
}
 
void SDPSolverInternal::evaluate(int nfdir, int nadir){
  throw CasadiException("SDPSolverInternal::evaluate: Not implemented");
}
 
void SDPSolverInternal::solve(){
  throw CasadiException("SDPSolverInternal::solve: Not implemented");
}
 
} // namespace CasADi

  


