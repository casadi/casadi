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


#include "csparse_cholesky_internal.hpp"

/// \cond INTERNAL

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINSOL_CSPARSECHOLESKY_EXPORT
  casadi_register_linsol_csparsecholesky(Linsol::Plugin* plugin) {
    plugin->creator = CSparseCholeskyInternal::creator;
    plugin->name = "csparsecholesky";
    plugin->doc = CSparseCholeskyInternal::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_CSPARSECHOLESKY_EXPORT casadi_load_linsol_csparsecholesky() {
    Linsol::registerPlugin(casadi_register_linsol_csparsecholesky);
  }

  CSparseCholeskyInternal::CSparseCholeskyInternal(const std::string& name,
                                                   const Sparsity& sparsity, int nrhs) :
    Linsol(name, sparsity, nrhs) {
    L_ = 0;
    S_ = 0;

    casadi_assert_message(sparsity.is_symmetric(),
                          "CSparseCholeskyInternal: supplied sparsity must be symmetric, got "
                          << sparsity.dim() << ".");
  }

  CSparseCholeskyInternal::~CSparseCholeskyInternal() {
    if (S_) cs_sfree(S_);
    if (L_) cs_nfree(L_);
  }

  void CSparseCholeskyInternal::init() {
    // Call the init method of the base class
    Linsol::init();

    A_.nzmax = nnz_in(0);  // maximum number of entries
    A_.m = size1_in(0); // number of columns
    A_.n = size2_in(0); // number of rows
    A_.p = const_cast<int*>(sparsity_in(0).colind()); // row pointers (size n+1)
                                                         // or row indices (size nzmax)
    A_.i = const_cast<int*>(sparsity_in(0).row()); // column indices, size nzmax
    A_.x = 0; // numerical values, size nzmax
    A_.nz = -1; // of entries in triplet matrix, -1 for compressed-row

    // Temporary
    temp_.resize(A_.n);

    if (verbose()) {
      userOut() << "CSparseCholeskyInternal::prepare: symbolic factorization" << endl;
    }

    // ordering and symbolic analysis
    int order = 0; // ordering?
    if (S_) cs_sfree(S_);
    S_ = cs_schol(order, &A_) ;
  }


  Sparsity CSparseCholeskyInternal::linsol_cholesky_sparsity(bool tr) const {
    casadi_assert(S_);
    int n = A_.n;
    int nzmax = S_->cp[n];
    std::vector< int > row(n+1);
    std::copy(S_->cp, S_->cp+n+1, row.begin());
    std::vector< int > colind(nzmax);
    int *Li = &colind.front();
    int *Lp = &row.front();
    const cs* C;
    C = S_->pinv ? cs_symperm(&A_, S_->pinv, 1) : &A_;
    std::vector< int > temp(2*n);
    int *c = & temp.front();
    int *s = c+n;
    for (int k = 0 ; k < n ; k++) c[k] = S_->cp[k] ;
    for (int k = 0 ; k < n ; k++) {       /* compute L(k, :) for L*L' = C */
      int top = cs_ereach(C, k, S_->parent, s, c) ;
      for ( ; top < n ; top++) {  /* solve L(0:k-1, 0:k-1) * x = C(:, k) */
          int i = s[top] ;               /* s[top..n-1] is pattern of L(k, :) */
          int p = c[i]++ ;
          Li[p] = k ;                /* store L(k, i) in row i */
      }
      int p = c[k]++ ;
      Li[p] = k ;
    }
    Lp[n] = S_->cp[n] ;
    Sparsity ret(n, n, row, colind); // BUG?

    return tr ? ret.T() : ret;

  }

  DM CSparseCholeskyInternal::linsol_cholesky(bool tr) const {
    casadi_assert(L_);
    cs *L = L_->L;
    int nz = L->nzmax;
    int m = L->m; // number of cols
    int n = L->n; // number of rows
    std::vector< int > colind(m+1);
    std::copy(L->p, L->p+m+1, colind.begin());
    std::vector< int > row(nz);
    std::copy(L->i, L->i+nz, row.begin());
    std::vector< double > data(nz);
    std::copy(L->x, L->x+nz, data.begin());
    DM ret(Sparsity(n, m, colind, row), data, false);

    return tr ? ret.T() : ret;
  }

  void CSparseCholeskyInternal::linsol_factorize(Memory& m, const double* A) {
    // Set the nonzeros of the matrix
    casadi_assert(A!=0);
    A_.x = const_cast<double*>(A);

    // Make sure that all entries of the linear system are valid
    int nnz = nnz_in(0);
    for (int k=0; k<nnz; ++k) {
      casadi_assert_message(!isnan(A[k]), "Nonzero " << k << " is not-a-number");
      casadi_assert_message(!isinf(A[k]), "Nonzero " << k << " is infinite");
    }

    if (L_) cs_nfree(L_);
    L_ = cs_chol(&A_, S_) ;                 // numeric Cholesky factorization
    casadi_assert(L_!=0);
  }

  void CSparseCholeskyInternal::linsol_solve(Memory& m, double* x, int nrhs, bool tr) {
    casadi_assert(L_!=0);

    double *t = &temp_.front();
    for (int k=0; k<nrhs; ++k) {
      if (tr) {
        cs_pvec(S_->q, x, t, A_.n) ;   // t = P1\b
        cs_ltsolve(L_->L, t) ;               // t = L\t
        cs_lsolve(L_->L, t) ;              // t = U\t
        cs_pvec(L_->pinv, t, x, A_.n) ;      // x = P2\t
      } else {
        cs_ipvec(L_->pinv, x, t, A_.n) ;   // t = P1\b
        cs_lsolve(L_->L, t) ;               // t = L\t
        cs_ltsolve(L_->L, t) ;              // t = U\t
        cs_ipvec(S_->q, t, x, A_.n) ;      // x = P2\t
      }
      x += ncol();
    }
  }

  void CSparseCholeskyInternal::linsol_solveL(double* x, int nrhs, bool tr) {
    casadi_assert(L_!=0);

    double *t = getPtr(temp_);

    for (int k=0; k<nrhs; ++k) {
      cs_ipvec(L_->pinv, x, t, A_.n) ;   // t = P1\b
      if (tr) cs_lsolve(L_->L, t) ; // t = L\t
      if (!tr) cs_ltsolve(L_->L, t) ; // t = U\t
      cs_ipvec(S_->q, t, x, A_.n) ;      // x = P2\t
      x += ncol();
    }
  }

} // namespace casadi

/// \endcond
