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

  void LapackQr::init() {
    // Call the base class initializer
    Linsol::init();

    // Get dimensions
    ncol_ = ncol();
    nrow_ = nrow();

    // Currently only square matrices tested
    if (ncol_!=nrow_) throw CasadiException("LapackQr::init: currently only "
                                           "square matrices implemented.");

    // Allocate matrix
    mat_.resize(ncol_*ncol_);
    tau_.resize(ncol_);
    work_.resize(10*ncol_);
  }

  void LapackQr::linsol_factorize(Memory& mem, const double* A) const {
    LapackQr& m = const_cast<LapackQr&>(*this);

    // Get the elements of the matrix, dense format
    casadi_densify(A, sparsity_, getPtr(m.mat_), false);

    // Factorize the matrix
    int info = -100;
    int lwork = work_.size();
    dgeqrf_(&m.ncol_, &m.ncol_, getPtr(m.mat_), &m.ncol_, getPtr(m.tau_),
            getPtr(m.work_), &lwork, &info);
    if (info != 0) throw CasadiException("LapackQr::prepare: dgeqrf_ "
                                         "failed to factorize the Jacobian");
  }

  void LapackQr::linsol_solve(Memory& mem, double* x, int nrhs, bool tr) const {
    LapackQr& m = const_cast<LapackQr&>(*this);

    // Properties of R
    char uploR = 'U';
    char diagR = 'N';
    char sideR = 'L';
    double alphaR = 1.;
    char transR = tr ? 'T' : 'N';

    // Properties of Q
    char transQ = tr ? 'N' : 'T';
    char sideQ = 'L';
    int k = tau_.size(); // minimum of ncol_ and nrow_
    int lwork = work_.size();

    if (tr) {

      // Solve for transpose(R)
      dtrsm_(&sideR, &uploR, &transR, &diagR, &m.ncol_, &nrhs, &alphaR,
             getPtr(m.mat_), &m.ncol_, x, &m.ncol_);

      // Multiply by Q
      int info = 100;
      dormqr_(&sideQ, &transQ, &m.ncol_, &nrhs, &k, getPtr(m.mat_), &m.ncol_, getPtr(m.tau_), x,
              &m.ncol_, getPtr(m.work_), &lwork, &info);
      if (info != 0) throw CasadiException("LapackQr::solve: dormqr_ failed "
                                          "to solve the linear system");

    } else {

      // Multiply by transpose(Q)
      int info = 100;
      dormqr_(&sideQ, &transQ, &m.ncol_, &nrhs, &k, getPtr(m.mat_), &m.ncol_, getPtr(m.tau_), x,
              &m.ncol_, getPtr(m.work_), &lwork, &info);
      if (info != 0) throw CasadiException("LapackQr::solve: dormqr_ failed to "
                                          "solve the linear system");

      // Solve for R
      dtrsm_(&sideR, &uploR, &transR, &diagR, &m.ncol_, &nrhs, &alphaR,
             getPtr(m.mat_), &m.ncol_, x, &m.ncol_);
    }
  }

} // namespace casadi
