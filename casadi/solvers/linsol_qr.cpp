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

  LinsolQr::LinsolQr(const std::string& name, const Sparsity& sp)
    : LinsolInternal(name, sp) {
  }

  LinsolQr::~LinsolQr() {
    clear_mem();
  }

  void LinsolQr::init(const Dict& opts) {
    // Call the init method of the base class
    LinsolInternal::init(opts);

    // Dimensions
    int nrow = this->nrow(), ncol = this->ncol();

    // Allocate memory
    leftmost_.resize(nrow);
    parent_.resize(ncol);
    pinv_.resize(nrow + ncol);
    vector<int> iw(nrow + 7*ncol + 1);

    // Initialize QP solve
    casadi_qr_init(sp_, sp_.T(),
                   get_ptr(leftmost_), get_ptr(parent_), get_ptr(pinv_),
                   &nrow_ext_, &v_nnz_, &r_nnz_, get_ptr(iw));

    // Calculate sparsities
    sp_v_.resize(2 + ncol + 1 + v_nnz_);
    sp_r_.resize(2 + ncol + 1 + r_nnz_);
    casadi_qr_sparsities(sp_, nrow_ext_, get_ptr(sp_v_), get_ptr(sp_r_),
                         get_ptr(leftmost_), get_ptr(parent_), get_ptr(pinv_),
                         get_ptr(iw));
  }

  int LinsolQr::init_mem(void* mem) const {
    if (LinsolInternal::init_mem(mem)) return 1;
    auto m = static_cast<LinsolQrMemory*>(mem);

    // Memory for numerical solution
    m->nz_v.resize(v_nnz_);
    m->nz_r.resize(r_nnz_);
    m->beta.resize(ncol());
    m->w.resize(nrow_ext_);
    m->iw.resize(nrow() + 7*ncol() + 1);
    m->y.resize(nrow_ext_);
    return 0;
  }

  int LinsolQr::sfact(void* mem, const double* A) const {
    return 0;
  }

  int LinsolQr::nfact(void* mem, const double* A) const {
    auto m = static_cast<LinsolQrMemory*>(mem);
    casadi_qr(sp_, A, get_ptr(m->iw), get_ptr(m->w),
              get_ptr(sp_v_), get_ptr(m->nz_v), get_ptr(sp_r_), get_ptr(m->nz_r),
              get_ptr(m->beta), get_ptr(leftmost_), get_ptr(parent_), get_ptr(pinv_));
    return 0;
  }

  int LinsolQr::solve(void* mem, const double* A, double* x, int nrhs, bool tr) const {
    auto m = static_cast<LinsolQrMemory*>(mem);
    int nrow_ext = m->y.size();
    int ncol = this->ncol();
    for (int k=0; k<nrhs; ++k) {
      if (tr) {

        // ('P'Q R)' x = R'Q'P x = b <-> x = P' Q R' \ b
        // Copy to y
        casadi_copy(x, ncol, get_ptr(m->y));
        //  Solve for R'
        casadi_qr_trs(get_ptr(sp_r_), get_ptr(m->nz_r), get_ptr(m->y), 1);
        // Multiply by Q
        casadi_qr_mv(get_ptr(sp_v_), get_ptr(m->nz_v), get_ptr(m->beta),
                     get_ptr(m->y), 0);
        // Multiply by P'
        for (int c=0; c<ncol; ++c) x[c] = m->y.at(pinv_.at(c));
      } else {
        //P'Q R x = b <-> x = R \ Q' P b
        // Multiply with P
        casadi_fill(get_ptr(m->y), nrow_ext, 0.);
        for (int c=0; c<ncol; ++c) m->y.at(pinv_.at(c)) = x[c];
        // Multiply with Q'
        casadi_qr_mv(get_ptr(sp_v_), get_ptr(m->nz_v), get_ptr(m->beta),
                     get_ptr(m->y), 1);
        //  Solve for R
        casadi_qr_trs(get_ptr(sp_r_), get_ptr(m->nz_r), get_ptr(m->y), 0);
        // Copy to x
        casadi_copy(get_ptr(m->y), ncol, x);
      }
      x += ncol;
    }

    return 0;
  }

  void LinsolQr::generate(CodeGenerator& g, const std::string& A, const std::string& x,
                          int nrhs, bool tr) const {
    // Not ready
    return LinsolInternal::generate(g, A, x, nrhs, tr);
    #if 0
    // Place in block to avoid conflicts caused by local variables
    g << "{";
    // Work vectors TODO(@jaeandersson): Use work vectors from Solve
    g << "real_t nz_v[" << nz_v.size() << "], "
         "nz_v[" << nz_v.size() << "], "
         "beta[" << beta.size() << "], "
         "w[" << w.size() << "], "
         "y[" << y.size() << "];\n";
    g << "int iw[" << iw.size() << "];\n";

    // End of block
    g << "}\n"
    #endif
  }

} // namespace casadi
