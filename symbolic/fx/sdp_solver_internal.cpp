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

// Constructor
SDPSolverInternal::SDPSolverInternal(const std::vector<CRSSparsity> &st) : st_(st) {
  addOption("calc_p",OT_BOOLEAN, true, "Indicate if the P-part of primal solution should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (m x m).");
  addOption("calc_dual",OT_BOOLEAN, true, "Indicate if dual should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (m x m).");

  casadi_assert_message(st_.size()==SDP_STRUCT_NUM,"Problem structure mismatch");
  
  const CRSSparsity& A = st_[SDP_STRUCT_A];
  const CRSSparsity& G = st_[SDP_STRUCT_G];
  const CRSSparsity& F = st_[SDP_STRUCT_F];
  
  casadi_assert_message(G==G.transpose(),"SDPSolverInternal: Supplied G sparsity must symmetric but got " << G.dimString());
  
  m_ = G.size1();
  
  nc_ = A.size1();
  n_ = A.size2();
  
  casadi_assert_message(F.size2()==m_,"SDPSolverInternal: Supplied F sparsity: number of columns (" << F.size2() <<  ")  must match m (" << m_ << ")");
  
  casadi_assert_message(F.size1()%n_==0,"SDPSolverInternal: Supplied F sparsity: number of rows (" << F.size2() <<  ")  must be an integer multiple of n (" << n_ << "), but got remainder " << F.size1()%n_);
  
  // Input arguments
  setNumInputs(SDP_SOLVER_NUM_IN);
  input(SDP_SOLVER_G) = DMatrix(G,0);
  input(SDP_SOLVER_F) = DMatrix(F,0);
  input(SDP_SOLVER_A) = DMatrix(A,0);
  input(SDP_SOLVER_C) = DMatrix::zeros(n_);
  input(SDP_SOLVER_LBX) = -DMatrix::inf(n_);
  input(SDP_SOLVER_UBX) = DMatrix::inf(n_);
  input(SDP_SOLVER_LBA) = -DMatrix::inf(nc_);
  input(SDP_SOLVER_UBA) = DMatrix::inf(nc_);

  for (int i=0;i<n_;i++) {
    CRSSparsity s = input(SDP_SOLVER_F)(range(i*m_,(i+1)*m_),ALL).sparsity();
    casadi_assert_message(s==s.transpose(),"SDPSolverInternal: Each supplied Fi must be symmetric. But got " << s.dimString() <<  " for i = " << i << ".");
  }
  
  inputScheme_ = SCHEME_SDPInput;
  outputScheme_ = SCHEME_SDPOutput;

}
    
void SDPSolverInternal::init() {
  // Call the init method of the base class
  FXInternal::init();
  
  calc_p_ = getOption("calc_p");
  calc_dual_ = getOption("calc_dual");

  // Find aggregate sparsity pattern
  CRSSparsity aggregate = input(SDP_SOLVER_G).sparsity();
  for (int i=0;i<n_;++i) {
    aggregate = aggregate + input(SDP_SOLVER_F)(range(i*m_,(i+1)*m_),ALL).sparsity();
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
  
  if (nb_>0) {
    // Make a mapping function from (G,F) -> (G[p,p]_j,F_i[p,p]j)
    SXMatrix G = ssym("G",input(SDP_SOLVER_G).sparsity());
    SXMatrix F = ssym("F",input(SDP_SOLVER_F).sparsity());

    std::vector<SXMatrix> in;
    in.push_back(G);
    in.push_back(F);
    std::vector<SXMatrix> out((n_+1)*nb_);
    for (int j=0;j<nb_;++j) {
      out[j] = G(p,p)(range(r[j],r[j+1]),range(r[j],r[j+1]));
    }
    for (int i=0;i<n_;++i) {
      SXMatrix Fi = F(range(i*m_,(i+1)*m_),ALL)(p,p);
      for (int j=0;j<nb_;++j) {
        out[(i+1)*nb_+j] = Fi(range(r[j],r[j+1]),range(r[j],r[j+1]));
      }
    }
    mapping_ = SXFunction(in,out);
    mapping_.init();
  }

  // Output arguments
  setNumOutputs(SDP_SOLVER_NUM_OUT);
  output(SDP_SOLVER_X) = DMatrix::zeros(n_,1);
  output(SDP_SOLVER_P) = calc_p_? DMatrix(Pmapper_.output().sparsity(),0) : DMatrix();
  output(SDP_SOLVER_DUAL) = calc_dual_? DMatrix(Pmapper_.output().sparsity(),0) : DMatrix();
  output(SDP_SOLVER_COST) = 0.0;
  output(SDP_SOLVER_DUAL_COST) = 0.0;
  output(SDP_SOLVER_LAM_X) = DMatrix::zeros(n_,1);
  output(SDP_SOLVER_LAM_A) = DMatrix::zeros(nc_,1);
  
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

  


