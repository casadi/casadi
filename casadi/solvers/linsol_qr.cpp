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
    int nrow_ext, v_nnz, r_nnz;
    casadi_qr_init(sp_, sp_.T(),
                   get_ptr(leftmost_), get_ptr(parent_), get_ptr(pinv_),
                   &nrow_ext, &v_nnz, &r_nnz, get_ptr(iw));

    // Calculate sparsities
    vector<int> sp_v(2 + ncol + 1 + v_nnz);
    vector<int> sp_r(2 + ncol + 1 + r_nnz);
    casadi_qr_sparsities(sp_, nrow_ext, get_ptr(sp_v), get_ptr(sp_r),
                         get_ptr(leftmost_), get_ptr(parent_), get_ptr(pinv_),
                         get_ptr(iw));
    sp_v_ = Sparsity::compressed(sp_v);
    sp_r_ = Sparsity::compressed(sp_r);
  }

  int LinsolQr::init_mem(void* mem) const {
    if (LinsolInternal::init_mem(mem)) return 1;
    auto m = static_cast<LinsolQrMemory*>(mem);

    // Memory for numerical solution
    m->v.resize(sp_v_.nnz());
    m->r.resize(sp_r_.nnz());
    m->beta.resize(ncol());
    m->w.resize(sp_r_.size1());
    m->iw.resize(nrow() + 7*ncol() + 1);
    return 0;
  }

  int LinsolQr::sfact(void* mem, const double* A) const {
    return 0;
  }

  int LinsolQr::nfact(void* mem, const double* A) const {
    auto m = static_cast<LinsolQrMemory*>(mem);
    casadi_qr(sp_, A, get_ptr(m->iw), get_ptr(m->w),
              sp_v_, get_ptr(m->v), sp_r_, get_ptr(m->r),
              get_ptr(m->beta), get_ptr(leftmost_), get_ptr(parent_), get_ptr(pinv_));
    return 0;
  }

  int LinsolQr::solve(void* mem, const double* A, double* x, int nrhs, bool tr) const {
    auto m = static_cast<LinsolQrMemory*>(mem);
    casadi_qr_solve(x, nrhs, tr,
                    sp_v_, get_ptr(m->v), sp_r_, get_ptr(m->r),
                    get_ptr(m->beta), get_ptr(pinv_), get_ptr(m->w));
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
    g << "real_t v[" << sp_v_.nnz() << "], "
         "r[" << sp_r_.nnz() << "], "
         "beta[" << ncol() << "], "
         "w[" << sp_r_.size1() << "], "
    g << "int iw[" << nrow() + 7*ncol() + 1 << "];\n";

    // End of block
    g << "}\n"
    #endif
  }

} // namespace casadi
