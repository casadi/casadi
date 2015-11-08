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


#include "csparse_interface.hpp"
#include "casadi/core/casadi_options.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINSOL_CSPARSE_EXPORT
  casadi_register_linsol_csparse(Linsol::Plugin* plugin) {
    plugin->creator = CsparseInterface::creator;
    plugin->name = "csparse";
    plugin->doc = CsparseInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_CSPARSE_EXPORT casadi_load_linsol_csparse() {
    Linsol::registerPlugin(casadi_register_linsol_csparse);
  }

  CsparseInterface::CsparseInterface(const std::string& name,
                                     const Sparsity& sparsity, int nrhs)
    : Linsol(name, sparsity, nrhs) {
    N_ = 0;
    S_ = 0;
  }

  CsparseInterface::~CsparseInterface() {
    if (S_) cs_sfree(S_);
    if (N_) cs_nfree(N_);
  }

  void CsparseInterface::init() {
    // Call the init method of the base class
    Linsol::init();

    A_.nzmax = input().nnz();  // maximum number of entries
    A_.m = input().size1(); // number of rows
    A_.n = input().size2(); // number of columns
    A_.p = const_cast<int*>(input().colind()); // column pointers (size n+1)
                                                        // or col indices (size nzmax)
    A_.i = const_cast<int*>(input().row()); // row indices, size nzmax
    A_.x = &input()->front(); // numerical values, size nzmax
    A_.nz = -1; // of entries in triplet matrix, -1 for compressed-col

    // Temporary
    temp_.resize(A_.n);

    // Has the routine been called once
    called_once_ = false;
  }

  void CsparseInterface::
  linsol_prepare(const double** arg, double** res, int* iw, double* w, void* mem) {
    Linsol::linsol_prepare(arg, res, iw, w, mem);

    if (!called_once_) {
      if (verbose()) {
        userOut() << "CsparseInterface::prepare: symbolic factorization" << endl;
      }

      // ordering and symbolic analysis
      int order = 0; // ordering?
      if (S_) cs_sfree(S_);
      S_ = cs_sqr(order, &A_, 0) ;
    }

    called_once_ = true;

    // Get a referebce to the nonzeros of the linear system
    const vector<double>& linsys_nz = input().data();

    // Make sure that all entries of the linear system are valid
    for (int k=0; k<linsys_nz.size(); ++k) {
      casadi_assert_message(!isnan(linsys_nz[k]), "Nonzero " << k << " is not-a-number");
      casadi_assert_message(!isinf(linsys_nz[k]), "Nonzero " << k << " is infinite");
    }

    if (verbose()) {
      userOut() << "CsparseInterface::prepare: numeric factorization" << endl;
      userOut() << "linear system to be factorized = " << endl;
      input(0).printSparse();
    }

    double tol = 1e-8;

    if (N_) cs_nfree(N_);
    N_ = cs_lu(&A_, S_, tol) ;                 // numeric LU factorization
    if (N_==0) {
      DMatrix temp = sparsify(input());
      if (temp.sparsity().is_singular()) {
        stringstream ss;
        ss << "CsparseInterface::prepare: factorization failed due to matrix"
          " being singular. Matrix contains numerical zeros which are "
            "structurally non-zero. Promoting these zeros to be structural "
            "zeros, the matrix was found to be structurally rank deficient."
            " sprank: " << sprank(temp.sparsity()) << " <-> " << temp.size2() << endl;
        if (verbose()) {
          ss << "Sparsity of the linear system: " << endl;
          input(LINSOL_A).sparsity().print(ss); // print detailed
        }
        throw CasadiException(ss.str());
      } else {
        stringstream ss;
        ss << "CsparseInterface::prepare: factorization failed, check if Jacobian is singular"
           << endl;
        if (verbose()) {
          ss << "Sparsity of the linear system: " << endl;
          input(LINSOL_A).sparsity().print(ss); // print detailed
        }
        throw CasadiException(ss.str());
      }
    }
    casadi_assert(N_!=0);
  }

  void CsparseInterface::linsol_solve(double* x, int nrhs, bool transpose) {
    double time_start=0;
    casadi_assert(N_!=0);

    double *t = &temp_.front();

    for (int k=0; k<nrhs; ++k) {
      if (transpose) {
        cs_pvec(S_->q, x, t, A_.n) ;       // t = P2*b
        casadi_assert(N_->U!=0);
        cs_utsolve(N_->U, t) ;              // t = U'\t
        cs_ltsolve(N_->L, t) ;              // t = L'\t
        cs_pvec(N_->pinv, t, x, A_.n) ;    // x = P1*t
      } else {
        cs_ipvec(N_->pinv, x, t, A_.n) ;   // t = P1\b
        cs_lsolve(N_->L, t) ;               // t = L\t
        cs_usolve(N_->U, t) ;               // t = U\t
        cs_ipvec(S_->q, t, x, A_.n) ;      // x = P2\t
      }
      x += ncol();
    }
  }

} // namespace casadi
