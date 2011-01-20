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

using namespace std;
namespace CasADi{

LapackQRDense::LapackQRDense(){
}

LapackQRDense::LapackQRDense(int nrow, int ncol, int nrhs){
  assignNode(new LapackQRDenseInternal(nrow,ncol,nrhs));
}
 
LapackQRDenseInternal* LapackQRDense::operator->(){
  return static_cast<LapackQRDenseInternal*>(FX::operator->());
}

const LapackQRDenseInternal* LapackQRDense::operator->() const{
  return static_cast<const LapackQRDenseInternal*>(FX::operator->());
}

LapackQRDenseInternal::LapackQRDenseInternal(int nrow, int ncol, int nrhs) : LinearSolverInternal(nrow,ncol,nrhs){
    
  // Currently only square matrices tested
  if(nrow!=ncol) throw CasadiException("LapackQRDenseInternal::LapackQRDenseInternal: currently only square matrices implemented.");
}

LapackQRDenseInternal::~LapackQRDenseInternal(){
}

void LapackQRDenseInternal::init(){
  // Call the base class initializer
  LinearSolverInternal::init();
  
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
  dgeqrf_(&nrow_, &nrow_, &mat_[0], &nrow_, &tau_[0], &work_[0], &lwork, &info);
  if(info != 0) throw CasadiException("LapackQRDenseInternal::prepare: dgeqrf_ failed to factorize the jacobian");

  // Success if reached this point
  prepared_ = true;
}
    
void LapackQRDenseInternal::solve(){
  // Input and output vectors
  const vector<double>& b = input(1);
  vector<double>& x = output();
  
  // Copy the right hand side to the solution vector
  copy(b.begin(),b.end(),x.begin());

  int info = 100;
  
  // Solve for transpose(R)
  char uplo = 'U';
  char transR = 'T';
  char diag = 'N';
  char sideR = 'L';
  double alpha = 1.;
  dtrsm_(&sideR, &uplo, &transR, &diag, &nrow_, &nrhs_, &alpha, &mat_[0], &nrow_, &x[0], &nrow_);
  
  // Multiply by Q
  char transQ = 'N';
  char sideQ = 'L';
  int k = tau_.size(); // minimum of nrow_ and ncol_
  int lwork = work_.size();
  dormqr_(&sideQ, &transQ, &nrow_, &nrhs_, &k, &mat_[0], &nrow_, &tau_[0], &x[0], &nrow_, &work_[0], &lwork, &info);
  if(info != 0) throw CasadiException("LapackQRDenseInternal::solve: dormqr_ failed to solve the linear system");
}


} // namespace CasADi

  


