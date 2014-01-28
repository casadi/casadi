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

#include "lapack_qr_nullspace_internal.hpp"
#include "../../symbolic/fx/fx_internal.hpp"


using namespace std;
namespace CasADi{

  
  LapackQRNullspaceInternal::LapackQRNullspaceInternal(const CRSSparsity& A_sp, int nfwd, int nadj) : NullspaceInternal(A_sp), nfwd_(nfwd), nadj_(nadj) {
 
    setOption("dense",true);
  
  }
  
  LapackQRNullspaceInternal::~LapackQRNullspaceInternal(){
  }


  void LapackQRNullspaceInternal::init(){
    NullspaceInternal::init();
    
    casadi_assert(dense_);
    setNumInputs((1+nfwd_) + nadj_);
    
    CRSSparsity dense = output(0).sparsity();
    
    for (int i=0;i<nfwd_+1;++i) {
      input(i)=DMatrix(A_sp_);
    }
    for (int i=0;i<nadj_;++i) {
      input((nfwd_+1)+i)=DMatrix(dense);
    }
    
    setNumOutputs((1+nfwd_) + nadj_);
    for (int i=0;i<nfwd_+1;++i) {
      output(i)=DMatrix(dense);
    }
    for (int i=0;i<nadj_;++i) {
      output((nfwd_+1)+i)=DMatrix(A_sp_);
    }
    
    tempA_.resize(n_*m_);
    
    tau_.resize(min(n_,m_));
    
    work_.resize(m_*64);
    work2_.resize(n_*64);
  
  }

  void LapackQRNullspaceInternal::evaluate(){
    // Copy input
    std::copy(input(0).data().begin(),input(0).data().end(),tempA_.begin());
    int info;
    
    int lwork = work_.size();
    
    dgeqrf_(&n_,&m_,&tempA_[0],&n_,&tau_[0],&work_[0],&lwork,&info);
    casadi_assert(info==0);
    
    char side = 'R';
    char trans = 'T';
    
    int nm = n_-m_;
    
    // Obtain the last nm columns of Q
    output(0).setAll(0.0);
    for (int i=0;i<nm;++i) {
      output(0).data()[nm*(i+m_)+i]=1;
    }

    int lwork2 = work2_.size();
    
    // 
    dormqr_(&side, &trans, &nm, &n_, &m_, &tempA_[0], &n_, &tau_[0], &output(0).data()[0], &nm, &work2_[0], &lwork2, &info);
    casadi_assert(info==0);

    std::cout << "out" << output(0) << std::endl;
  
  }
  
  LapackQRNullspaceInternal* LapackQRNullspaceInternal::clone() const{
    // Return a deep copy
    LapackQRNullspaceInternal* node = new LapackQRNullspaceInternal(A_sp_, nfwd_, nadj_);
    node->setOption(dictionary());
    return node;
  }
  
  FX LapackQRNullspaceInternal::getDerivative(int nfwd, int nadj) {
    casadi_assert(nfwd_==0 && nadj_==0);
    
    LapackQRNullspaceInternal* node = new LapackQRNullspaceInternal(A_sp_,nfwd, nadj);
    node->setOption(dictionary());
    
    LapackQRNullspace ret;
    ret.assignNode(node);

    return ret;
    
  }


} // namespace CasADi


