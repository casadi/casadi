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
#include "casadi/core/matrix/matrix_tools.hpp"

/// \cond INTERNAL

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINEARSOLVER_CSPARSECHOLESKY_EXPORT
  casadi_register_linearsolver_csparsecholesky(LinearSolverInternal::Plugin* plugin) {
    plugin->creator = CSparseCholeskyInternal::creator;
    plugin->name = "csparsecholesky";
    plugin->doc = CSparseCholeskyInternal::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_LINEARSOLVER_CSPARSECHOLESKY_EXPORT casadi_load_linearsolver_csparsecholesky() {
    LinearSolverInternal::registerPlugin(casadi_register_linearsolver_csparsecholesky);
  }

  CSparseCholeskyInternal::CSparseCholeskyInternal(const Sparsity& sparsity, int nrhs) :
      LinearSolverInternal(sparsity, nrhs) {
    L_ = 0;
    S_ = 0;

    casadi_assert_message(sparsity.isSymmetric(),
                          "CSparseCholeskyInternal: supplied sparsity must be symmetric, got "
                          << sparsity.dimString() << ".");
  }

  CSparseCholeskyInternal::CSparseCholeskyInternal(const CSparseCholeskyInternal& linsol) :
      LinearSolverInternal(linsol) {
    L_ = 0;
    S_ = 0;
  }

  CSparseCholeskyInternal::~CSparseCholeskyInternal() {
    if (S_) cs_sfree(S_);
    if (L_) cs_nfree(L_);
  }

  void CSparseCholeskyInternal::init() {
    // Call the init method of the base class
    LinearSolverInternal::init();

    AT_.nzmax = input().size();  // maximum number of entries
    AT_.m = input().size1(); // number of cols
    AT_.n = input().size2(); // number of rows
    AT_.p = const_cast<int*>(&input().colind().front()); // row pointers (size n+1)
                                                         // or row indices (size nzmax)
    AT_.i = const_cast<int*>(&input().row().front()); // col indices, size nzmax
    AT_.x = &input().front(); // col indices, size nzmax
    AT_.nz = -1; // of entries in triplet matrix, -1 for compressed-row

    // Temporary
    temp_.resize(AT_.n);

    if (verbose()) {
      cout << "CSparseCholeskyInternal::prepare: symbolic factorization" << endl;
    }

    // ordering and symbolic analysis
    int order = 0; // ordering?
    if (S_) cs_sfree(S_);
    S_ = cs_schol(order, &AT_) ;
  }


  Sparsity CSparseCholeskyInternal::getFactorizationSparsity(bool transpose) const {
    casadi_assert(S_);
    int n = AT_.n;
    int nzmax = S_->cp[n];
    std::vector< int > row(n+1);
    std::copy(S_->cp, S_->cp+n+1, row.begin());
    std::vector< int > colind(nzmax);
    int *Li = &colind.front();
    int *Lp = &row.front();
    const cs* C;
    C = S_->pinv ? cs_symperm(&AT_, S_->pinv, 1) : &AT_;
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

    return transpose? ret.T() : ret;

  }

  DMatrix CSparseCholeskyInternal::getFactorization(bool transpose) const {
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
    DMatrix ret(Sparsity(n, m, colind, row), data);

    return transpose? ret.T() : ret;
  }

  void CSparseCholeskyInternal::prepare() {

    prepared_ = false;

    // Get a reference to the nonzeros of the linear system
    const vector<double>& linsys_nz = input().data();

    // Make sure that all entries of the linear system are valid
    for (int k=0; k<linsys_nz.size(); ++k) {
      casadi_assert_message(!isnan(linsys_nz[k]), "Nonzero " << k << " is not-a-number");
      casadi_assert_message(!isinf(linsys_nz[k]), "Nonzero " << k << " is infinite");
    }

    if (verbose()) {
      cout << "CSparseCholeskyInternal::prepare: numeric factorization" << endl;
      cout << "linear system to be factorized = " << endl;
      input(0).printSparse();
    }

    if (L_) cs_nfree(L_);
    L_ = cs_chol(&AT_, S_) ;                 // numeric Cholesky factorization
    if (L_==0) {
      DMatrix temp = input();
      temp.sparsify();
      if (temp.sparsity().isSingular()) {
        stringstream ss;
        ss << "CSparseCholeskyInternal::prepare: factorization failed due "
          "to matrix being singular. Matrix contains numerical zeros which are"
          " structurally non-zero. Promoting these zeros to be structural "
          "zeros, the matrix was found to be structurally rank deficient. "
          "sprank: " << rank(temp.sparsity()) << " <-> " << temp.size2() << endl;
        if (verbose()) {
          ss << "Sparsity of the linear system: " << endl;
          input(LINSOL_A).sparsity().print(ss); // print detailed
        }
        throw CasadiException(ss.str());
      } else {
        stringstream ss;
        ss << "CSparseCholeskyInternal::prepare: factorization failed, "
            "check if Jacobian is singular" << endl;
        if (verbose()) {
          ss << "Sparsity of the linear system: " << endl;
          input(LINSOL_A).sparsity().print(ss); // print detailed
        }
        throw CasadiException(ss.str());
      }
    }
    casadi_assert(L_!=0);

    prepared_ = true;
  }

  void CSparseCholeskyInternal::solve(double* x, int nrhs, bool transpose) {
    casadi_assert(prepared_);
    casadi_assert(L_!=0);

    double *t = &temp_.front();
    for (int k=0; k<nrhs; ++k) {
      if (transpose) {
        cs_pvec(S_->q, x, t, AT_.n) ;   // t = P1\b
        cs_ltsolve(L_->L, t) ;               // t = L\t
        cs_lsolve(L_->L, t) ;              // t = U\t
        cs_pvec(L_->pinv, t, x, AT_.n) ;      // x = P2\t
      } else {
        cs_ipvec(L_->pinv, x, t, AT_.n) ;   // t = P1\b
        cs_lsolve(L_->L, t) ;               // t = L\t
        cs_ltsolve(L_->L, t) ;              // t = U\t
        cs_ipvec(S_->q, t, x, AT_.n) ;      // x = P2\t
      }
      x += ncol();
    }
  }

  void CSparseCholeskyInternal::solveL(double* x, int nrhs, bool transpose) {
    casadi_assert(prepared_);
    casadi_assert(L_!=0);

    double *t = getPtr(temp_);

    for (int k=0; k<nrhs; ++k) {
      cs_ipvec(L_->pinv, x, t, AT_.n) ;   // t = P1\b
      if (transpose) cs_lsolve(L_->L, t) ; // t = L\t
      if (!transpose) cs_ltsolve(L_->L, t) ; // t = U\t
      cs_ipvec(S_->q, t, x, AT_.n) ;      // x = P2\t
      x += ncol();
    }
  }

  CSparseCholeskyInternal* CSparseCholeskyInternal::clone() const {
    return new CSparseCholeskyInternal(input(LINSOL_A).sparsity(), input(LINSOL_B).size2());
  }

} // namespace casadi

/// \endcond
