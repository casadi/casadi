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
#include "../../symbolic/std_vector_tools.hpp"

using namespace std;
namespace CasADi{

  LapackQRDense::LapackQRDense(){
  }

  LapackQRDense::LapackQRDense(const Sparsity& sparsity, int nrhs){
    assignNode(new LapackQRDenseInternal(sparsity,nrhs));
  }
 
  LapackQRDenseInternal* LapackQRDense::operator->(){
    return static_cast<LapackQRDenseInternal*>(Function::operator->());
  }

  const LapackQRDenseInternal* LapackQRDense::operator->() const{
    return static_cast<const LapackQRDenseInternal*>(Function::operator->());
  }

  LapackQRDenseInternal::LapackQRDenseInternal(const Sparsity& sparsity, int nrhs) : LinearSolverInternal(sparsity,nrhs){
  }

  LapackQRDenseInternal::~LapackQRDenseInternal(){
  }

  void LapackQRDenseInternal::init(){
    // Call the base class initializer
    LinearSolverInternal::init();
  
    // Get dimensions
    ncol_ = ncol();
    nrow_ = nrow();
  
    // Currently only square matrices tested
    if(ncol_!=nrow_) throw CasadiException("LapackQRDenseInternal::init: currently only square matrices implemented.");
  
    // Allocate matrix
    mat_.resize(ncol_*ncol_);
    tau_.resize(ncol_);
    work_.resize(10*ncol_);
  }

  void LapackQRDenseInternal::prepare(){
    prepared_ = false;
  
    // Get the elements of the matrix, dense format
    input(0).get(mat_,DENSE);
  
    // Factorize the matrix
    int info = -100;
    int lwork = work_.size();
    dgeqrf_(&ncol_, &ncol_, getPtr(mat_), &ncol_, getPtr(tau_), getPtr(work_), &lwork, &info);
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
    int k = tau_.size(); // minimum of ncol_ and nrow_
    int lwork = work_.size();
  
    if(transpose){

      // Solve for transpose(R)
      dtrsm_(&sideR, &uploR, &transR, &diagR, &ncol_, &nrhs, &alphaR, getPtr(mat_), &ncol_, x, &ncol_);
    
      // Multiply by Q
      int info = 100;
      dormqr_(&sideQ, &transQ, &ncol_, &nrhs, &k, getPtr(mat_), &ncol_, getPtr(tau_), x, &ncol_, getPtr(work_), &lwork, &info);
      if(info != 0) throw CasadiException("LapackQRDenseInternal::solve: dormqr_ failed to solve the linear system");

    } else {

      // Multiply by transpose(Q)
      int info = 100;
      dormqr_(&sideQ, &transQ, &ncol_, &nrhs, &k, getPtr(mat_), &ncol_, getPtr(tau_), x, &ncol_, getPtr(work_), &lwork, &info);
      if(info != 0) throw CasadiException("LapackQRDenseInternal::solve: dormqr_ failed to solve the linear system");

      // Solve for R
      dtrsm_(&sideR, &uploR, &transR, &diagR, &ncol_, &nrhs, &alphaR, getPtr(mat_), &ncol_, x, &ncol_);  
    }
  }

  LapackQRDenseInternal* LapackQRDenseInternal::clone() const{
    return new LapackQRDenseInternal(*this);
  }

} // namespace CasADi
