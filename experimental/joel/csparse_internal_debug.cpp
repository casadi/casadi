/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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
#include "core/matrix/crs_sparsity_internal.hpp"
#include "core/std_vector_tools.hpp"

using namespace std;
namespace casadi{

CSparseInternal::CSparseInternal(const CRSSparsity& sparsity)  : LinearSolverInternal(sparsity){
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
    int order = 3; // ordering?
    int qr = 1; // LU
    if(S_) cs_sfree(S_);
    S_ = cs_sqr (order, &AT_, qr) ;              
    
    std::vector<int> pinv;
    std::vector<int> q;
    std::vector<int> parent; 
    std::vector<int> cp;
    std::vector<int> leftmost;
    int m2;
    double lnz;
    double unz;
    input().sparsity()->prefactorize(order, qr, pinv, q, parent, cp, leftmost, m2, lnz, unz);
    

    cout << "pinv" << endl;
    cout << pinv << endl;
    if(S_->pinv!=0)
      cout << vector<int>(S_->pinv, S_->pinv + pinv.size()) << endl;
    cout << endl;
      
    cout << "q" << endl;
    cout << q << endl;
    if(S_->q!=0)
      cout << vector<int>(S_->q, S_->q + q.size()) << endl;
    cout << endl;

    cout << "parent" << endl;
    cout << parent << endl;
    if(S_->parent!=0)
      cout << vector<int>(S_->parent, S_->parent + parent.size()) << endl;
    cout << endl;

    cout << "cp" << endl;
    cout << cp << endl;
    if(S_->cp!=0)
      cout << vector<int>(S_->cp, S_->cp + cp.size()) << endl;
    cout << endl;
    
    cout << "leftmost" << endl;
    cout << leftmost << endl;
    if(S_->leftmost!=0)
      cout << vector<int>(S_->leftmost, S_->leftmost + leftmost.size()) << endl;
    cout << endl;
    
    
    
  }
  
  prepared_ = false;
  called_once_ = true;

  double tol = 1e-8;
  
  if(N_) cs_nfree(N_);
  N_ = cs_lu(&AT_, S_, tol) ;                 // numeric LU factorization 
  if(N_==0){
    throw CasadiException("factorization failed, Jacobian singular?");
  }
  casadi_assert(N_!=0);

  prepared_ = true;
}
  
void CSparseInternal::solve(double* x, int nrhs, bool transpose){
  casadi_assert(prepared_);
  casadi_assert(N_!=0);
  
  double *t = &temp_.front();
  
  for(int k=0; k<nrhs; ++k){
    if(transpose){
      cs_ipvec (N_->pinv, x, t, AT_.n) ;   // t = P1\b
      cs_lsolve (N_->L, t) ;               // t = L\t 
      cs_usolve (N_->U, t) ;               // t = U\t 
      cs_ipvec (S_->q, t, x, AT_.n) ;      // x = P2\t 
    } else {
      cs_pvec (S_->q, x, t, AT_.n) ;       // t = P2*b 
      casadi_assert(N_->U!=0);
      cs_utsolve (N_->U, t) ;              // t = U'\t 
      cs_ltsolve (N_->L, t) ;              // t = L'\t 
      cs_pvec (N_->pinv, t, x, AT_.n) ;    // x = P1*t 
    }
    x += nrow();
  }
}


CSparseInternal* CSparseInternal::clone() const{
  return new CSparseInternal(sparsity_);
}

} // namespace casadi
