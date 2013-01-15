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

#include "dsdp_internal.hpp"

#include "../../symbolic/stl_vector_tools.hpp"
#include "../../symbolic/matrix/matrix_tools.hpp"

/**
Some implementation details

"Multiple cones can be created for the same solver, but it is usually more efficient to group all blocks into the same conic structure." user manual

*/
using namespace std;
namespace CasADi {

DSDPInternal* DSDPInternal::clone() const{
  // Return a deep copy
  DSDPInternal* node = new DSDPInternal(input(SDP_C).sparsity(),input(SDP_A).sparsity());
  if(!node->is_init_)
    node->init();
  return node;
}
  
DSDPInternal::DSDPInternal(const CRSSparsity &C, const CRSSparsity &A) : SDPSolverInternal(C,A){
 
  casadi_assert_message(double(n_)*(double(n_)+1)/2 < std::numeric_limits<int>::max(),"Your problem size n is too large to be handled by DSDP.");
    
}

DSDPInternal::~DSDPInternal(){ 

}

void DSDPInternal::init(){
  SDPSolverInternal::init();


  // A return flag used by DSDP
  int info;
  
  // Allocate DSDP solver memory
  info = DSDPCreate(m_, &dsdp);
  info = DSDPCreateSDPCone(dsdp,1,&sdpcone);
  info = SDPConeSetBlockSize(sdpcone, 0, n_);
  info = SDPConeSetSparsity(sdpcone, 0, n_);

  // Fill the data structures that hold DSDP-style sparse symmetric matrix
  pattern_.resize(m_+1);
  values_.resize(m_+1);
  
  pattern_[0].resize(input(SDP_C).sparsity().sizeL());
  values_[0].resize(pattern_[0].size());
  int nz=0;
  vector<int> rowind,col;
  input(SDP_C).sparsity().getSparsityCRS(rowind,col);
  for(int r=0; r<rowind.size()-1; ++r) {
    for(int el=rowind[r]; el<rowind[r+1]; ++el){
     if(r>=col[el]){
       pattern_[0][nz++] = r*(r + 1)/2 + col[el];
     }
    }
  }
  
  for (int i=0;i<m_;++i) {
    CRSSparsity Ai = mapping_.output(i).sparsity();
    pattern_[i+1].resize(Ai.sizeL());
    values_[i+1].resize(pattern_[i+1].size());
    int nz=0;
    vector<int> rowind,col;
    Ai.getSparsityCRS(rowind,col);
    for(int r=0; r<rowind.size()-1; ++r) {
      for(int el=rowind[r]; el<rowind[r+1]; ++el){
       if(r>=col[el]){
         pattern_[i+1][nz++] = r*(r + 1)/2 + col[el];
       }
      }
    }
    mapping_.output(i).get(values_[i+1],SPARSESYM);
  }
  
  if (calc_dual_) {
    store_X_.resize(n_*(n_+1)/2);
  }
  if (calc_p_) {
    store_P_.resize(n_*(n_+1)/2);
  }
}

void DSDPInternal::evaluate(int nfdir, int nadir) {
  int info;
  
  // Copy b vector
  for (int i=0;i<m_;++i) {
    info = DSDPSetDualObjective(dsdp, i+1, -input(SDP_B).at(i));
  }
  
  // Get Ai from supplied A
  mapping_.setInput(input(SDP_A));
  // Negate because the standard form in PSDP is different
  std::transform(mapping_.input().begin(), mapping_.input().end(), mapping_.input().begin(), std::negate<double>());
  mapping_.evaluate();

  // Set 
  input(SDP_C).get(values_[0],SPARSESYM);
  std::transform(values_[0].begin(), values_[0].end(), values_[0].begin(), std::negate<double>());
  info = SDPConeSetASparseVecMat(sdpcone, 0, 0, n_, 1, 0, &pattern_[0][0], &values_[0][0], pattern_[0].size() );
  
  for (int i=0;i<m_;++i) {
    mapping_.output(i).get(values_[i+1],SPARSESYM);
    info = SDPConeSetASparseVecMat(sdpcone, 0, i+1, n_, 1, 0, &pattern_[i+1][0], &values_[i+1][0], pattern_[i+1].size() );
  }
  
  info = DSDPSetup(dsdp);
  
  info = DSDPSolve(dsdp);
  casadi_assert_message(info==0,"DSDPSolver failed");
  
  info = DSDPGetY(dsdp,&output(SDP_PRIMAL).at(0),m_);
  
  double temp;
  DSDPGetDDObjective(dsdp, &temp);
  output(SDP_PRIMAL_COST).set(-temp);
  DSDPGetPPObjective(dsdp, &temp);
  output(SDP_DUAL_COST).set(-temp);
  
  if (calc_dual_) {
    info = SDPConeComputeX(sdpcone, 0, n_, &store_X_[0], store_X_.size());
    output(SDP_DUAL).set(store_X_,SPARSESYM);
  }
  //

  if (calc_p_) {
    info = SDPConeComputeS(sdpcone, 0, 1.0,  &output(SDP_PRIMAL).at(0), m_, 0, n_ , &store_P_[0], store_P_.size());
    output(SDP_PRIMAL_P).set(store_P_,SPARSESYM);
  }
  

  
}

} // namespace CasADi
