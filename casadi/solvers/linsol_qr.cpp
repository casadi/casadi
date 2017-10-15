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


#include "linsol_qr.hpp"
#include "casadi/core/global_options.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINSOL_QR_EXPORT
  casadi_register_linsol_qr(LinsolInternal::Plugin* plugin) {
    plugin->creator = LinsolQr::creator;
    plugin->name = "qr";
    plugin->doc = LinsolQr::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &LinsolQr::options_;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_QR_EXPORT casadi_load_linsol_qr() {
    LinsolInternal::registerPlugin(casadi_register_linsol_qr);
  }

  LinsolQr::LinsolQr(const std::string& name)
    : LinsolInternal(name) {
  }

  LinsolQr::~LinsolQr() {
    clear_mem();
  }

  LinsolQrMemory::~LinsolQrMemory() {
  }

  void LinsolQr::init(const Dict& opts) {
    // Call the init method of the base class
    LinsolInternal::init(opts);
  }

  int LinsolQr::init_mem(void* mem) const {
    return LinsolInternal::init_mem(mem);
  }

  void LinsolQr::reset(void* mem, const int* sp) const {
    LinsolInternal::reset(mem, sp);
    auto m = static_cast<LinsolQrMemory*>(mem);


    // Temporary
//    m->temp_.resize(m->A.n);
  }

  void LinsolQr::pivoting(void* mem, const double* A) const {
    LinsolInternal::pivoting(mem, A);
    auto m = static_cast<LinsolQrMemory*>(mem);

   }

  void LinsolQr::factorize(void* mem, const double* A) const {
    auto m = static_cast<LinsolQrMemory*>(mem);

#if 0
    // Set the nonzeros of the matrix
    m->A.x = const_cast<double*>(A);

{
    auto S2 = cs_sqr(0, &m->A, 1);
    auto N2 = cs_qr(&m->A, S2);

    // Get V
    int V_nrow = N2->L->m;
    int V_ncol = N2->L->n;
    const int *V_colind = N2->L->p;
    const int *V_row = N2->L->i;
    Sparsity V_sp(V_nrow, V_ncol, V_colind, V_row);
    const double* V_nz = N2->L->x;
    DM V(V_sp, vector<double>(V_nz, V_nz+V_sp.nnz()));
    uout() << "V = " << V << endl;

    // Get R
    int R_nrow = N2->U->m;
    int R_ncol = N2->U->n;
    const int *R_colind = N2->U->p;
    const int *R_row = N2->U->i;
    Sparsity R_sp(R_nrow, R_ncol, R_colind, R_row);
    const double* R_nz = N2->U->x;
    DM R(R_sp, vector<double>(R_nz, R_nz+R_sp.nnz()));
    uout() << "R = " << R << endl;


//S2->pinv


    Sparsity Asp = Sparsity::compressed(m->sparsity);
    int size1=m->nrow(), size2=m->ncol();

    // Dimensions
    cout << "A.dim() = " << Asp.dim() << endl;

    // Allocate vectors
    std::vector<int> iw;
    std::vector<int> parent(size2);
    std::vector<int> post(size2);
    std::vector<int> count(size2);

    // Allocate work
    iw.resize(size1+size2);
    // Calculate elimination tree
    casadi_etree(Asp, get_ptr(parent), get_ptr(iw), true);
    // Calculate postorder
    iw.resize(3*size2);
    casadi_postorder(get_ptr(parent), size2, get_ptr(post), get_ptr(iw));
    // Calculate colind in L
    iw.resize(size1 + 5*size2 + 1);
    std::vector<int> L_colind(1+size2);
    casadi_qr_colind(Asp.T(), get_ptr(parent), get_ptr(post),
                     get_ptr(L_colind), get_ptr(iw));


    // Number of nonzeros in R:
    int r_nnz = L_colind.back();
    cout << "r_nnz = " << r_nnz << ", S2->unz = " << S2->unz << endl;

    // Number of nonzeros in V:
    iw.resize(size1 + 3*size2);
    vector<int> pinv(size1 + size2);
    vector<int> leftmost(size1);
    int nrow_ext = -111;
    cout << "parent = " << parent << endl;
    cout << "S2->parent " << vector<int>(S2->parent, S2->parent+size2) << endl;


    int v_nnz = casadi_qr_nnz(Asp, get_ptr(pinv), get_ptr(leftmost),
                              get_ptr(parent), &nrow_ext, get_ptr(iw));
    cout << "nrow_ext = " << nrow_ext << ", S2->m2 = " << S2->m2 << endl;

    cout << "v_nnz = " << v_nnz << ", S2->lnz = " <<  S2->lnz << endl;
    cout << "pinv = " << pinv << endl;
    cout << "S2->pinv = " << vector<int>(S2->pinv, S2->pinv+size1+size2) << endl;
    cout << "leftmost = " << leftmost << endl;
    cout << "S2->leftmost = " << vector<int>(S2->leftmost, S2->leftmost+size1) << endl;


    vector<double> nz_v(v_nnz, -1), nz_r(r_nnz, -1), beta(size2);
    vector<int> sp_v(2 + size2 + 1 + v_nnz, -11);
    vector<int> sp_r(2 + size2 + 1 + r_nnz, -22);
    sp_v[0] = sp_r[0] = nrow_ext;
    sp_v[1] = sp_r[1] = size2;
    iw.resize(nrow_ext + size2);
    vector<double> w(nrow_ext);
    casadi_qr(Asp, A, get_ptr(iw), get_ptr(w),
                   get_ptr(sp_v), get_ptr(nz_v), get_ptr(sp_r), get_ptr(nz_r),
                   get_ptr(beta), get_ptr(leftmost), get_ptr(parent), get_ptr(pinv));

    DM V2(Sparsity::compressed(sp_v), nz_v);
    uout() << "V2 = " << V2 << endl;
    DM R2(Sparsity::compressed(sp_r), nz_r);
    uout() << "R2 = " << R2 << endl;


    cout << "sp_v " << sp_v << endl;
    cout << "sp_r " << sp_r << endl;

    for (int c=0; c<size2; ++c) {
      vector<double> q1(nrow_ext, 0);
      q1.at(pinv.at(c)) = 1;
      casadi_qr_mv(get_ptr(sp_v), get_ptr(nz_v), get_ptr(beta),
                   get_ptr(q1), 1);

//      vector<double> q2(size1, -99);
//      for (int k=0; k<size1; ++k) {
  //      q2.at(k) = q1.at(pinv.at(k));
  //    }

      cout << "q2 = " << q1 << endl;
    }


    DM A3(Asp, vector<double>(A, A+Asp.nnz()));
    DM Q3, R3;
    qr(A3, Q3, R3);
    cout << "Q3 = " << Q3 << endl;
    cout << "R3 = " << R3 << endl;


    cs_nfree(N2);
    cs_sfree(S2);
}

    // Make sure that all entries of the linear system are valid
    for (int k=0; k<m->nnz(); ++k) {
      casadi_assert(!isnan(A[k]),
        "Nonzero " + str(k) + " is not-a-number");
      casadi_assert(!isinf(A[k]),
        "Nonzero " + str(k) + " is infinite");
    }

    if (verbose_) {
      uout() << "LinsolQr::prepare: numeric factorization" << endl;
      uout() << "linear system to be factorized = " << endl;
      Sparsity sp = Sparsity::compressed(m->sparsity);
      DM(sp, vector<double>(A, A+m->nnz())).print_sparse(uout());
    }

    double tol = 1e-8;

    if (m->N) cs_nfree(m->N);
    m->N = cs_lu(&m->A, m->S, tol) ;                 // numeric LU factorization
    if (m->N==0) {
      Sparsity sp = Sparsity::compressed(m->sparsity);
      DM temp(sp, vector<double>(A, A+sp.nnz()));
      temp = sparsify(temp);
      if (temp.sparsity().is_singular()) {
        stringstream ss;
        ss << "LinsolQr::prepare: factorization failed due to matrix"
          " being singular. Matrix contains numerical zeros which are "
            "structurally non-zero. Promoting these zeros to be structural "
            "zeros, the matrix was found to be structurally rank deficient."
            " sprank: " << sprank(temp.sparsity()) << " <-> " << temp.size2() << endl;
        if (verbose_) {
          ss << "Sparsity of the linear system: " << endl;
          sp.disp(ss, true); // print detailed
        }
        throw CasadiException(ss.str());
      } else {
        stringstream ss;
        ss << "LinsolQr::prepare: factorization failed, check if Jacobian is singular"
           << endl;
        if (verbose_) {
          ss << "Sparsity of the linear system: " << endl;
          sp.disp(ss, true); // print detailed
        }
        throw CasadiException(ss.str());
      }
    }
    casadi_assert_dev(m->N!=0);
#endif
  }

  void LinsolQr::solve(void* mem, double* x, int nrhs, bool tr) const {
    auto m = static_cast<LinsolQrMemory*>(mem);
    //casadi_assert_dev(m->N!=0);

#if 0
    double *t = &m->temp_.front();

    for (int k=0; k<nrhs; ++k) {
      if (tr) {
        cs_pvec(m->S->q, x, t, m->A.n) ;       // t = P2*b
        casadi_assert_dev(m->N->U!=0);
        cs_utsolve(m->N->U, t) ;              // t = U'\t
        cs_ltsolve(m->N->L, t) ;              // t = L'\t
        cs_pvec(m->N->pinv, t, x, m->A.n) ;    // x = P1*t
      } else {
        cs_ipvec(m->N->pinv, x, t, m->A.n) ;   // t = P1\b
        cs_lsolve(m->N->L, t) ;               // t = L\t
        cs_usolve(m->N->U, t) ;               // t = U\t
        cs_ipvec(m->S->q, t, x, m->A.n) ;      // x = P2\t
      }
      x += m->ncol();
    }
    #endif
  }

} // namespace casadi
