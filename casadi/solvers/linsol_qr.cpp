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

    // Dimensions
    int size1 = this->nrow(), size2 = this->ncol();

    // Allocate memory
    m->leftmost.resize(size1);
    m->parent.resize(size2);
    m->pinv.resize(size1 + size2);
    m->iw.resize(size1 + 7*size2 + 1);

    // Initialize QP solve
    int nrow_ext, v_nnz, r_nnz;
    casadi_qr_init(sp_, sp_.T(),
                   get_ptr(m->leftmost), get_ptr(m->parent), get_ptr(m->pinv),
                   &nrow_ext, &v_nnz, &r_nnz, get_ptr(m->iw));

    // Calculate sparsities
    m->sp_v.resize(2 + size2 + 1 + v_nnz);
    m->sp_r.resize(2 + size2 + 1 + r_nnz);
    casadi_qr_sparsities(sp_, nrow_ext, get_ptr(m->sp_v), get_ptr(m->sp_r),
                         get_ptr(m->leftmost), get_ptr(m->parent), get_ptr(m->pinv),
                         get_ptr(m->iw));

    // Memory for numerical solution
    m->nz_v.resize(v_nnz);
    m->nz_r.resize(r_nnz);
    m->beta.resize(size2);
    m->w.resize(nrow_ext);
    m->y.resize(nrow_ext);
  }

  void LinsolQr::pivoting(void* mem, const double* A) const {
    LinsolInternal::pivoting(mem, A);
    //auto m = static_cast<LinsolQrMemory*>(mem);
   }

  void LinsolQr::factorize(void* mem, const double* A) const {
    auto m = static_cast<LinsolQrMemory*>(mem);
    casadi_qr(sp_, A, get_ptr(m->iw), get_ptr(m->w),
              get_ptr(m->sp_v), get_ptr(m->nz_v), get_ptr(m->sp_r), get_ptr(m->nz_r),
              get_ptr(m->beta), get_ptr(m->leftmost), get_ptr(m->parent), get_ptr(m->pinv));
  }

  void LinsolQr::solve(void* mem, double* x, int nrhs, bool tr) const {
    auto m = static_cast<LinsolQrMemory*>(mem);
    int nrow_ext = m->y.size();
    int ncol = this->ncol();

    for (int k=0; k<nrhs; ++k) {
      if (tr) {
        // ('P'Q R)' x = R'Q'P x = b <-> x = P' Q R' \ b
        // Copy to y
        casadi_copy(x, ncol, get_ptr(m->y));
        //  Solve for R'
        casadi_qr_trs(get_ptr(m->sp_r), get_ptr(m->nz_r), get_ptr(m->y), 1);
        // Multiply by Q
        casadi_qr_mv(get_ptr(m->sp_v), get_ptr(m->nz_v), get_ptr(m->beta),
                     get_ptr(m->y), 0);
        // Multiply by P'
        for (int c=0; c<ncol; ++c) x[c] = m->y.at(m->pinv.at(c));
      } else {
        //P'Q R x = b <-> x = R \ Q' P b
        // Multiply with P
        casadi_fill(get_ptr(m->y), nrow_ext, 0.);
        for (int c=0; c<ncol; ++c) m->y.at(m->pinv.at(c)) = x[c];
        // Multiply with Q'
        casadi_qr_mv(get_ptr(m->sp_v), get_ptr(m->nz_v), get_ptr(m->beta),
                     get_ptr(m->y), 1);
        //  Solve for R
        casadi_qr_trs(get_ptr(m->sp_r), get_ptr(m->nz_r), get_ptr(m->y), 0);
        // Copy to x
        casadi_copy(get_ptr(m->y), ncol, x);
      }
      x += ncol;
    }
  }

} // namespace casadi
