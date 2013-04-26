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
#include "symbolic/matrix/matrix_tools.hpp"

using namespace std;
namespace CasADi{

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
        if(verbose()){
          cout << "CSparseInternal::prepare: symbolic factorization" << endl;
        }
        
        // ordering and symbolic analysis 
    int order = 0; // ordering?
    if(S_) cs_sfree(S_);
    S_ = cs_sqr (order, &AT_, 0) ;              
  }
  
  prepared_ = false;
  called_once_ = true;
  
  // Get a referebce to the nonzeros of the linear system
  const vector<double>& linsys_nz = input().data();
  
  // Make sure that all entries of the linear system are valid
  for(int k=0; k<linsys_nz.size(); ++k){
        casadi_assert_message(!isnan(linsys_nz[k]),"Nonzero " << k << " is not-a-number");
        casadi_assert_message(!isinf(linsys_nz[k]),"Nonzero " << k << " is infinite");
  }
  
  if(verbose()){
        cout << "CSparseInternal::prepare: numeric factorization" << endl;
        cout << "linear system to be factorized = " << endl;
        input(0).printSparse();
  }

  double tol = 1e-8;
  
  if(N_) cs_nfree(N_);
  N_ = cs_lu(&AT_, S_, tol) ;                 // numeric LU factorization 
  if(N_==0){
    DMatrix temp = input();
    makeSparse(temp);
    if (isSingular(temp.sparsity())) {
            stringstream ss;
            ss << "CSparseInternal::prepare: factorization failed due to matrix being singular. Matrix contains numerical zeros which are structurally non-zero. Promoting these zeros to be structural zeros, the matrix was found to be structurally rank deficient. sprank: " << rank(temp.sparsity()) << " <-> " << temp.size1() << endl;
            if(verbose()){
              ss << "Sparsity of the linear system: " << endl;
              sparsity_.print(ss); // print detailed
            }
      throw CasadiException(ss.str());
    } else {
            stringstream ss;
            ss << "CSparseInternal::prepare: factorization failed, check if Jacobian is singular" << endl;
            if(verbose()){
              ss << "Sparsity of the linear system: " << endl;
              sparsity_.print(ss); // print detailed
            }
      throw CasadiException(ss.str());
    }
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

} // namespace CasADi
