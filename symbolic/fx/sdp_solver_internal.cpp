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
  
  inputScheme_ = SCHEME_SDPInput;
  outputScheme_ = SCHEME_SDPOutput;
}



// Constructor
SDPSolverInternal::SDPSolverInternal(const CRSSparsity &C, const CRSSparsity &A) {
  addOption("calc_p",OT_BOOLEAN, true, "Indicate if the P-part of primal solution should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (n x n).");
  addOption("calc_dual",OT_BOOLEAN, true, "Indicate if dual should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (n x n).");
  
  casadi_assert_message(C==C.transpose(),"SDPSolverInternal: Supplied C sparsity must symmetric but got " << C.dimString());
  
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

  // Find aggregate sparsity pattern
  CRSSparsity aggregate = input(SDP_C).sparsity();
  for (int i=0;i<m_;++i) {
    aggregate = aggregate + input(SDP_A)(range(i*n_,(i+1)*n_),ALL).sparsity();
  }
  
  // Detect block diagonal structure in this sparsity pattern
  std::vector<int> p;
  std::vector<int> r;
  nb_ = aggregate.stronglyConnectedComponents(p,r);
  block_boundaries_.resize(nb_+1);
  std::copy(r.begin(),r.begin()+nb_+1,block_boundaries_.begin());  
  
  block_sizes_.resize(nb_);
  for (int i=0;i<nb_;++i) {
    block_sizes_[i]=r[i+1]-r[i];
  }
  
  // Make a mapping function from dense blocks to inversely-permuted block diagonal P
  std::vector< SXMatrix > full_blocks;
  for (int i=0;i<nb_;++i) {
    full_blocks.push_back(ssym("block",block_sizes_[i],block_sizes_[i]));
  }
  
  Pmapper_ = SXFunction(full_blocks,blkdiag(full_blocks)(lookupvector(p,p.size()),lookupvector(p,p.size())));
  Pmapper_.init();
  
  // Make a mapping function from (C,A) -> (C[p,p]_j,A_i[p,p]j)
  SXMatrix C = ssym("C",input(SDP_C).sparsity());
  SXMatrix A = ssym("A",input(SDP_A).sparsity());

  std::vector<SXMatrix> in;
  in.push_back(C);
  in.push_back(A);
  std::vector<SXMatrix> out((m_+1)*nb_);
  for (int j=0;j<nb_;++j) {
    out[j] = C(p,p)(range(r[j],r[j+1]),range(r[j],r[j+1]));
  }
  for (int i=0;i<m_;++i) {
    SXMatrix Ai = A(range(i*n_,(i+1)*n_),ALL)(p,p);
    for (int j=0;j<nb_;++j) {
      out[(i+1)*nb_+j] = Ai(range(r[j],r[j+1]),range(r[j],r[j+1]));
    }
  }
  mapping_ = SXFunction(in,out);
  mapping_.init();

  // Output arguments
  setNumOutputs(SDP_NUM_OUT);
  output(SDP_PRIMAL) = DMatrix::zeros(m_,1);
  output(SDP_PRIMAL_P) = calc_p_? DMatrix(Pmapper_.output().sparsity(),0) : DMatrix();
  output(SDP_DUAL) = calc_dual_? DMatrix(Pmapper_.output().sparsity(),0) : DMatrix();
  output(SDP_PRIMAL_COST) = 0.0;
  output(SDP_DUAL_COST) = 0.0;
  
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

  


