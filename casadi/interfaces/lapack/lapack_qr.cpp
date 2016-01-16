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
  casadi_register_linsol_lapackqr(Linsol::Plugin* plugin) {
    plugin->creator = LapackQr::creator;
    plugin->name = "lapackqr";
    plugin->doc = LapackQr::meta_doc.c_str();;
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_LAPACKQR_EXPORT casadi_load_linsol_lapackqr() {
    Linsol::registerPlugin(casadi_register_linsol_lapackqr);
  }

  LapackQr::LapackQr(const std::string& name,
                               const Sparsity& sparsity, int nrhs) :
    Linsol(name, sparsity, nrhs) {
  }

  LapackQr::~LapackQr() {
  }

  void LapackQr::init(const Dict& opts) {
    // Call the base class initializer
    Linsol::init(opts);

    // Currently only square matrices tested
    if (ncol()!=nrow()) throw CasadiException("LapackQr::init: currently only "
                                              "square matrices implemented.");
  }

  Memory* LapackQr::memory() const {
    LapackQrMemory* m = new LapackQrMemory();
    try {
      m->mat.resize(ncol()*ncol());
      m->tau.resize(ncol());
      m->work.resize(10*ncol());
      return m;
    } catch (...) {
      delete m;
      return 0;
    }
  }

  void LapackQr::linsol_factorize(Memory& mem, const double* A) const {
    LapackQrMemory& m = dynamic_cast<LapackQrMemory&>(mem);

    // Dimensions
    //int nrow = this->nrow();
    int ncol = this->ncol();

    // Get the elements of the matrix, dense format
    casadi_densify(A, sparsity_, get_ptr(m.mat), false);

    // Factorize the matrix
    int info = -100;
    int lwork = m.work.size();
    dgeqrf_(&ncol, &ncol, get_ptr(m.mat), &ncol, get_ptr(m.tau),
            get_ptr(m.work), &lwork, &info);
    if (info != 0) throw CasadiException("LapackQr::prepare: dgeqrf_ "
                                         "failed to factorize the Jacobian");
  }

  void LapackQr::linsol_solve(Memory& mem, double* x, int nrhs, bool tr) const {
    LapackQrMemory& m = dynamic_cast<LapackQrMemory&>(mem);

    // Dimensions
    int nrow = this->nrow();
    int ncol = this->ncol();

    // Properties of R
    char uploR = 'U';
    char diagR = 'N';
    char sideR = 'L';
    double alphaR = 1.;
    char transR = tr ? 'T' : 'N';

    // Properties of Q
    char transQ = tr ? 'N' : 'T';
    char sideQ = 'L';
    int k = m.tau.size(); // minimum of ncol and nrow
    int lwork = m.work.size();

    if (tr) {

      // Solve for transpose(R)
      dtrsm_(&sideR, &uploR, &transR, &diagR, &ncol, &nrhs, &alphaR,
             get_ptr(m.mat), &ncol, x, &ncol);

      // Multiply by Q
      int info = 100;
      dormqr_(&sideQ, &transQ, &ncol, &nrhs, &k, get_ptr(m.mat), &ncol, get_ptr(m.tau), x,
              &ncol, get_ptr(m.work), &lwork, &info);
      if (info != 0) throw CasadiException("LapackQr::solve: dormqr_ failed "
                                          "to solve the linear system");

    } else {

      // Multiply by transpose(Q)
      int info = 100;
      dormqr_(&sideQ, &transQ, &ncol, &nrhs, &k, get_ptr(m.mat), &ncol, get_ptr(m.tau), x,
              &ncol, get_ptr(m.work), &lwork, &info);
      if (info != 0) throw CasadiException("LapackQr::solve: dormqr_ failed to "
                                          "solve the linear system");

      // Solve for R
      dtrsm_(&sideR, &uploR, &transR, &diagR, &ncol, &nrhs, &alphaR,
             get_ptr(m.mat), &ncol, x, &ncol);
    }
  }

} // namespace casadi
