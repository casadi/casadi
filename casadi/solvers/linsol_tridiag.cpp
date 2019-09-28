/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *    Copyright (C) 2019 Jorn Baayen, KISTERS AG
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

#include <iostream>
#include "linsol_tridiag.hpp"
#include "casadi/core/global_options.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINSOL_TRIDIAG_EXPORT
  casadi_register_linsol_tridiag(LinsolInternal::Plugin* plugin) {
    plugin->creator = LinsolTridiag::creator;
    plugin->name = "tridiag";
    plugin->doc = LinsolTridiag::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &LinsolTridiag::options_;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_TRIDIAG_EXPORT casadi_load_linsol_tridiag() {
    LinsolInternal::registerPlugin(casadi_register_linsol_tridiag);
  }

  LinsolTridiag::LinsolTridiag(const std::string& name, const Sparsity& sp)
    : LinsolInternal(name, sp) {
  }

  LinsolTridiag::~LinsolTridiag() {
    clear_mem();
  }

  void LinsolTridiag::init(const Dict& opts) {
    // Call the init method of the base class
    LinsolInternal::init(opts);
  }

  int LinsolTridiag::init_mem(void* mem) const {
    if (LinsolInternal::init_mem(mem)) return 1;
    auto m = static_cast<LinsolTridiagMemory*>(mem);

    // Memory for numerical solution
    m->c.resize(nrow());
    m->ctr.resize(nrow());
    m->d.resize(nrow());
    return 0;
  }

  int LinsolTridiag::sfact(void* mem, const double* A) const {
    return 0;
  }

  int LinsolTridiag::nfact(void* mem, const double* A) const {
    auto m = static_cast<LinsolTridiagMemory*>(mem);
    m->have_c = false;
    m->have_ctr = false;
    return 0;
  }

  int LinsolTridiag::solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const {
    auto m = static_cast<LinsolTridiagMemory*>(mem);
    casadi_int i, k;
    const casadi_int *a_colind;
    // Extract sparsity
    a_colind = sp_ + 2;
    // Precompute cacheable coefficients
    if (tr) {
      if (!m->have_ctr) {
        // Precompute coefficients (transposed matrix)
        m->ctr[0] = A[a_colind[0] + 1] / A[a_colind[0] + 0];
        for (i = 1; i < nrow(); i++) {
          // A[i + 1, i] / (A[i, i] - A[i - 1, i] * m->ctr[i - 1]);
          double denom = A[a_colind[i] + 1] - A[a_colind[i] + 0] * m->ctr[i - 1];
          m->ctr[i] = A[a_colind[i] + 2] / denom;
        }
        m->have_ctr = true;
      }
    } else {
      if (!m->have_c) {
        // Precompute coefficients
        m->c[0] = A[a_colind[1] + 0] / A[a_colind[0] + 0];
        double denom = A[a_colind[1] + 1] - A[a_colind[1 - 1] + 1] * m->c[1 - 1];
        m->c[1] = A[a_colind[1 + 1] + 0] / denom;
        for (i = 2; i < nrow(); i++) {
          // A[i, i + 1] / (A[i, i] - A[i, i - 1] * m->c[i - 1]);
          double denom = A[a_colind[i] + 1] - A[a_colind[i - 1] + 2] * m->c[i - 1];
          m->c[i] = A[a_colind[i + 1] + 0] / denom;
        }
        m->have_c = true;
      }
    }
    // Compute
    for (k = 0; k < nrhs; k++) {
      if (tr) {
        // Compute remaining coefficients (transposed matrix)
        m->d[0] = x[0] / A[a_colind[0] + 0];
        for (i = 1; i < nrow(); i++) {
          //m->d[i] = (x[i] - A[i - 1, i] * m->d[i - 1]) / (A[i, i] - A[i - 1, i] * m->c[i - 1]);
          double denom = A[a_colind[i] + 1] - A[a_colind[i] + 0] * m->ctr[i - 1];
          m->d[i] = (x[i] - A[a_colind[i] + 0] * m->d[i - 1]) / denom;
        }
        // Solve
        x[nrow() - 1] = m->d[nrow() - 1];
        for (i = nrow() - 2; i >= 0; i--) {
          x[i] = m->d[i] - m->ctr[i] * x[i + 1];
        }
      } else {
        // Compute remaining coefficients
        m->d[0] = x[0] / A[a_colind[0] + 0];
        double denom = A[a_colind[1] + 1] - A[a_colind[1 - 1] + 1] * m->c[1 - 1];
        m->d[1] = (x[1] - A[a_colind[1 - 1] + 1] * m->d[1 - 1]) / denom;
        for (i = 2; i < nrow(); i++) {
          //m->d[i] = (x[i] - A[i, i - 1] * m->d[i - 1]) / (A[i, i] - A[i, i - 1] * m->c[i - 1]);
          double denom = A[a_colind[i] + 1] - A[a_colind[i - 1] + 2] * m->c[i - 1];
          m->d[i] = (x[i] - A[a_colind[i - 1] + 2] * m->d[i - 1]) / denom;
        }
        // Solve
        x[nrow() - 1] = m->d[nrow() - 1];
        for (i = nrow() - 2; i >= 0; i--) {
          x[i] = m->d[i] - m->c[i] * x[i + 1];
        }
      }
      x += ncol();
    }
    return 0;
  }

  void LinsolTridiag::generate(CodeGenerator& g, const std::string& A, const std::string& x,
                          casadi_int nrhs, bool tr) const {
    casadi_error("Not implemented");
  }

} // namespace casadi
