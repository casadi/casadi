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

#include "lapack_qr_dense.hpp"
#include "../../symbolic/stl_vector_tools.hpp"

using namespace std;
namespace CasADi{

  LapackQRDense::LapackQRDense(){
  }

  LapackQRDense::LapackQRDense(const CRSSparsity& sparsity, int nrhs){
    assignNode(new LapackQRDenseInternal(sparsity,nrhs));
  }
 
  LapackQRDenseInternal* LapackQRDense::operator->(){
    return static_cast<LapackQRDenseInternal*>(FX::operator->());
  }

  const LapackQRDenseInternal* LapackQRDense::operator->() const{
    return static_cast<const LapackQRDenseInternal*>(FX::operator->());
  }

  LapackQRDenseInternal::LapackQRDenseInternal(const CRSSparsity& sparsity, int nrhs) : LinearSolverInternal(sparsity,nrhs){
  }

  LapackQRDenseInternal::~LapackQRDenseInternal(){
  }

  void LapackQRDenseInternal::init(){
    // Call the base class initializer
    LinearSolverInternal::init();
  
    // Get dimensions
    nrow_ = nrow();
    ncol_ = ncol();
  
    // Currently only square matrices tested
    if(nrow_!=ncol_) throw CasadiException("LapackQRDenseInternal::init: currently only square matrices implemented.");
  
    // Allocate matrix
    mat_.resize(nrow_*nrow_);
    tau_.resize(nrow_);
    work_.resize(10*nrow_);
  }

  void LapackQRDenseInternal::prepare(){
    prepared_ = false;
  
    // Get the elements of the matrix, dense format
    input(0).get(mat_,DENSE);
  
    // Factorize the matrix
    int info = -100;
    int lwork = work_.size();
    dgeqrf_(&nrow_, &nrow_, getPtr(mat_), &nrow_, getPtr(tau_), getPtr(work_), &lwork, &info);
    if(info != 0) throw CasadiException("LapackQRDenseInternal::prepare: dgeqrf_ failed to factorize the jacobian");

    // Success if reached this point
    prepared_ = true;
  }
    
  void LapackQRDenseInternal::solve(double* x, int nrhs, bool transpose){
    // Properties of R
    char uploR = 'U';
    char diagR = 'N';
    char sideR = 'L';
    double alphaR = 1.;
    char transR = transpose ? 'T' : 'N';
  
    // Properties of Q
    char transQ = transpose ? 'N' : 'T';
    char sideQ = 'L';
    int k = tau_.size(); // minimum of nrow_ and ncol_
    int lwork = work_.size();
  
    if(transpose){

      // Solve for transpose(R)
      dtrsm_(&sideR, &uploR, &transR, &diagR, &nrow_, &nrhs, &alphaR, getPtr(mat_), &nrow_, x, &nrow_);
    
      // Multiply by Q
      int info = 100;
      dormqr_(&sideQ, &transQ, &nrow_, &nrhs, &k, getPtr(mat_), &nrow_, getPtr(tau_), x, &nrow_, getPtr(work_), &lwork, &info);
      if(info != 0) throw CasadiException("LapackQRDenseInternal::solve: dormqr_ failed to solve the linear system");

    } else {

      // Multiply by transpose(Q)
      int info = 100;
      dormqr_(&sideQ, &transQ, &nrow_, &nrhs, &k, getPtr(mat_), &nrow_, getPtr(tau_), x, &nrow_, getPtr(work_), &lwork, &info);
      if(info != 0) throw CasadiException("LapackQRDenseInternal::solve: dormqr_ failed to solve the linear system");

      // Solve for R
      dtrsm_(&sideR, &uploR, &transR, &diagR, &nrow_, &nrhs, &alphaR, getPtr(mat_), &nrow_, x, &nrow_);  
    }
  }

  LapackQRDenseInternal* LapackQRDenseInternal::clone() const{
    return new LapackQRDenseInternal(*this);
  }

} // namespace CasADi
