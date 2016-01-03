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
  }

  CsparseInterface::~CsparseInterface() {
  }

  void CsparseInterface::init() {
    // Call the init method of the base class
    Linsol::init();
  }

  void* CsparseInterface::alloc_mem() {
    Memory* m = new Memory();
    m->N = 0;
    m->S = 0;
    m->A.nzmax = nnz_in(0);  // maximum number of entries
    m->A.m = size1_in(0); // number of rows
    m->A.n = size2_in(0); // number of columns
    m->A.p = const_cast<int*>(sparsity_in(0).colind()); // column pointers (size n+1)
                                                     // or column indices (size nzmax)
    m->A.i = const_cast<int*>(sparsity_in(0).row()); // row indices, size nzmax
    m->A.x = 0; // numerical values, size nzmax
    m->A.nz = -1; // of entries in triplet matrix, -1 for compressed-col

    // Temporary
    m->temp_.resize(m->A.n);

    // Has the routine been called once
    m->called_once_ = false;

    return m;
  }

  void CsparseInterface::linsol_factorize(const double* A, void* mem) const {
    Memory* m = static_cast<Memory*>(mem);
    casadi_assert(A!=0);

    // Set the nonzeros of the matrix
    m->A.x = const_cast<double*>(A);

    if (!m->called_once_) {
      if (verbose()) {
        userOut() << "CsparseInterface::prepare: symbolic factorization" << endl;
      }

      // ordering and symbolic analysis
      int order = 0; // ordering?
      if (m->S) cs_sfree(m->S);
      m->S = cs_sqr(order, &m->A, 0) ;
    }

    m->called_once_ = true;

    // Make sure that all entries of the linear system are valid
    for (int k=0; k<sparsity_.nnz(); ++k) {
      casadi_assert_message(!isnan(A[k]), "Nonzero " << k << " is not-a-number");
      casadi_assert_message(!isinf(A[k]), "Nonzero " << k << " is infinite");
    }

    if (verbose()) {
      userOut() << "CsparseInterface::prepare: numeric factorization" << endl;
      userOut() << "linear system to be factorized = " << endl;
      DM(sparsity_, vector<double>(A, A+sparsity_.nnz())).printSparse();
    }

    double tol = 1e-8;

    if (m->N) cs_nfree(m->N);
    m->N = cs_lu(&m->A, m->S, tol) ;                 // numeric LU factorization
    if (m->N==0) {
      DM temp(sparsity_, vector<double>(A, A+sparsity_.nnz()));
      temp = sparsify(temp);
      if (temp.sparsity().is_singular()) {
        stringstream ss;
        ss << "CsparseInterface::prepare: factorization failed due to matrix"
          " being singular. Matrix contains numerical zeros which are "
            "structurally non-zero. Promoting these zeros to be structural "
            "zeros, the matrix was found to be structurally rank deficient."
            " sprank: " << sprank(temp.sparsity()) << " <-> " << temp.size2() << endl;
        if (verbose()) {
          ss << "Sparsity of the linear system: " << endl;
          sparsity_.print(ss); // print detailed
        }
        throw CasadiException(ss.str());
      } else {
        stringstream ss;
        ss << "CsparseInterface::prepare: factorization failed, check if Jacobian is singular"
           << endl;
        if (verbose()) {
          ss << "Sparsity of the linear system: " << endl;
          sparsity_.print(ss); // print detailed
        }
        throw CasadiException(ss.str());
      }
    }
    casadi_assert(m->N!=0);
  }

  void CsparseInterface::linsol_solve(double* x, int nrhs, bool tr, void* mem) const {
    Memory* m = static_cast<Memory*>(mem);
    casadi_assert(m->N!=0);

    double *t = &m->temp_.front();

    for (int k=0; k<nrhs; ++k) {
      if (tr) {
        cs_pvec(m->S->q, x, t, m->A.n) ;       // t = P2*b
        casadi_assert(m->N->U!=0);
        cs_utsolve(m->N->U, t) ;              // t = U'\t
        cs_ltsolve(m->N->L, t) ;              // t = L'\t
        cs_pvec(m->N->pinv, t, x, m->A.n) ;    // x = P1*t
      } else {
        cs_ipvec(m->N->pinv, x, t, m->A.n) ;   // t = P1\b
        cs_lsolve(m->N->L, t) ;               // t = L\t
        cs_usolve(m->N->U, t) ;               // t = U\t
        cs_ipvec(m->S->q, t, x, m->A.n) ;      // x = P2\t
      }
      x += ncol();
    }
  }

  void CsparseInterface::free_mem(void* mem) {
    Memory* m = static_cast<Memory*>(mem);
    if (m->S) cs_sfree(m->S);
    if (m->N) cs_nfree(m->N);
    delete m;
  }

} // namespace casadi
