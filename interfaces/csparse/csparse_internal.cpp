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

#include "csparse_internal.hpp"

using namespace std;
namespace CasADi{
  namespace Interfaces{

CSparseInternal::CSparseInternal(const CRSSparsity& sparsity, int nrhs)  : LinearSolverInternal(sparsity,nrhs){
  N_ = 0;
  S_ = 0;
}

CSparseInternal::CSparseInternal(const CSparseInternal& linsol) : LinearSolverInternal(linsol){
}

CSparseInternal::~CSparseInternal(){
  if(S_) cs_sfree(S_);
  if(N_) cs_nfree(N_);
}

void CSparseInternal::init(){
  // Call the init method of the base class
  LinearSolverInternal::init();
  casadi_assert(nrhs_==1);

  AT_.nzmax = input().size();  // maximum number of entries 
  AT_.m = input().size2(); // number of rows
  AT_.n = input().size1(); // number of columns
  AT_.p = const_cast<int*>(&input().rowind().front()); // column pointers (size n+1) or col indices (size nzmax)
  AT_.i = const_cast<int*>(&input().col().front()); // row indices, size nzmax
  AT_.x = &input().front(); // row indices, size nzmax
  AT_.nz = -1; // of entries in triplet matrix, -1 for compressed-col 

  // Temporary
  temp_.resize(AT_.n);
  
  // Has the routine been called once
  called_once_ = false;
}

void CSparseInternal::prepare(){
  if(!called_once_){
    // ordering and symbolic analysis 
    int order = 0; // ordering?
    if(S_) cs_sfree(S_);
    S_ = cs_sqr (order, &AT_, 0) ;              
  }
  
  prepared_ = false;
  called_once_ = true;

  double tol = 1e-8;
  
  if(N_) cs_nfree(N_);
  N_ = cs_lu (&AT_, S_, tol) ;                 // numeric LU factorization 

  prepared_ = true;
}
  
void CSparseInternal::solve(){
  
  const double *b = &input(1).front();
  double *t = &temp_.front();
  double *x = &output().front();
  
  if(transpose_){
    cs_ipvec (N_->pinv, b, t, AT_.n) ;   // t = P1\b
    cs_lsolve (N_->L, t) ;               // t = L\t 
    cs_usolve (N_->U, t) ;               // t = U\t 
    cs_ipvec (S_->q, t, x, AT_.n) ;      // x = P2\t 
  } else {
    cs_pvec (S_->q, b, t, AT_.n) ;       // t = P2*b 
    cs_utsolve (N_->U, t) ;              // t = U'\t 
    cs_ltsolve (N_->L, t) ;              // t = L'\t 
    cs_pvec (N_->pinv, t, x, AT_.n) ;    // x = P1*t 
  }
}


CSparseInternal* CSparseInternal::clone() const{
  return new CSparseInternal(sparsity_,nrhs_);
}

  } // namespace Interfaces
} // namespace CasADi
