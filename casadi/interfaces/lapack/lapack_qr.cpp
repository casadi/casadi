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


#include "lapack_qr.hpp"
#include "../../core/std_vector_tools.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINSOL_LAPACKQR_EXPORT
  casadi_register_linsol_lapackqr(LinsolInternal::Plugin* plugin) {
    plugin->creator = LapackQr::creator;
    plugin->name = "lapackqr";
    plugin->doc = LapackQr::meta_doc.c_str();;
    plugin->version = 30;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_LAPACKQR_EXPORT casadi_load_linsol_lapackqr() {
    LinsolInternal::registerPlugin(casadi_register_linsol_lapackqr);
  }

  LapackQr::LapackQr(const std::string& name, const Sparsity& sparsity) :
    LinsolInternal(name, sparsity) {
  }

  LapackQr::~LapackQr() {
    clear_memory();
  }

  void LapackQr::init(const Dict& opts) {
    // Call the base class initializer
    LinsolInternal::init(opts);
  }

  void LapackQr::init_memory(void* mem) const {
    LinsolInternal::init_memory(mem);
    if (!sparsity_.is_null()) reset(mem, sparsity_);
  }

  void LapackQr::reset(void* mem, const int* sp) const {
    LinsolInternal::reset(mem, sp);
    auto m = static_cast<LapackQrMemory*>(mem);
    m->mat.resize(m->ncol() * m->ncol());
    m->tau.resize(m->ncol());
    m->work.resize(10*m->ncol());
  }

  void LapackQr::factorize(void* mem, const double* A) const {
    auto m = static_cast<LapackQrMemory*>(mem);

    // Dimensions
    //int nrow = this->nrow();
    int ncol = m->ncol();

    // Get the elements of the matrix, dense format
    casadi_densify(A, get_ptr(m->sparsity), get_ptr(m->mat), false);

    // Factorize the matrix
    int info = -100;
    int lwork = m->work.size();
    dgeqrf_(&ncol, &ncol, get_ptr(m->mat), &ncol, get_ptr(m->tau),
            get_ptr(m->work), &lwork, &info);
    if (info != 0) throw CasadiException("LapackQr::prepare: dgeqrf_ "
                                         "failed to factorize the Jacobian");
  }

  void LapackQr::solve(void* mem, double* x, int nrhs, bool tr) const {
    auto m = static_cast<LapackQrMemory*>(mem);

    // Dimensions
    //int nrow = this->nrow();
    int ncol = m->ncol();

    // Properties of R
    char uploR = 'U';
    char diagR = 'N';
    char sideR = 'L';
    double alphaR = 1.;
    char transR = tr ? 'T' : 'N';

    // Properties of Q
    char transQ = tr ? 'N' : 'T';
    char sideQ = 'L';
    int k = m->tau.size(); // minimum of ncol and nrow
    int lwork = m->work.size();

    if (tr) {

      // Solve for transpose(R)
      dtrsm_(&sideR, &uploR, &transR, &diagR, &ncol, &nrhs, &alphaR,
             get_ptr(m->mat), &ncol, x, &ncol);

      // Multiply by Q
      int info = 100;
      dormqr_(&sideQ, &transQ, &ncol, &nrhs, &k, get_ptr(m->mat), &ncol, get_ptr(m->tau), x,
              &ncol, get_ptr(m->work), &lwork, &info);
      if (info != 0) throw CasadiException("LapackQr::solve: dormqr_ failed "
                                          "to solve the linear system");

    } else {

      // Multiply by transpose(Q)
      int info = 100;
      dormqr_(&sideQ, &transQ, &ncol, &nrhs, &k, get_ptr(m->mat), &ncol, get_ptr(m->tau), x,
              &ncol, get_ptr(m->work), &lwork, &info);
      if (info != 0) throw CasadiException("LapackQr::solve: dormqr_ failed to "
                                          "solve the linear system");

      // Solve for R
      dtrsm_(&sideR, &uploR, &transR, &diagR, &ncol, &nrhs, &alphaR,
             get_ptr(m->mat), &ncol, x, &ncol);
    }
  }

} // namespace casadi
