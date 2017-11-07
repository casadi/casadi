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


#include "linsol_ldl.hpp"
#include "casadi/core/global_options.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINSOL_LDL_EXPORT
  casadi_register_linsol_ldl(LinsolInternal::Plugin* plugin) {
    plugin->creator = LinsolLdl::creator;
    plugin->name = "ldl";
    plugin->doc = LinsolLdl::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &LinsolLdl::options_;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_LDL_EXPORT casadi_load_linsol_ldl() {
    LinsolInternal::registerPlugin(casadi_register_linsol_ldl);
  }

  LinsolLdl::LinsolLdl(const std::string& name, const Sparsity& sp)
    : LinsolInternal(name, sp) {
  }

  LinsolLdl::~LinsolLdl() {
    clear_mem();
  }

  LinsolLdlMemory::~LinsolLdlMemory() {
  }

  void LinsolLdl::init(const Dict& opts) {
    // Call the init method of the base class
    LinsolInternal::init(opts);
  }

  int LinsolLdl::init_mem(void* mem) const {
    if (LinsolInternal::init_mem(mem)) return 1;
    auto m = static_cast<LinsolLdlMemory*>(mem);

    // Dimension
    int n=this->nrow();
    // Work vector
    m->iw.resize(3*n);
    // Elimination tree
    m->parent.resize(n);
    // Calculate colind in L (strictly lower entries only)
    m->sp_l.resize(2+n+1);
    m->sp_l[0] = n;
    m->sp_l[1] = n;
    casadi_ldl_colind(sp_, get_ptr(m->parent), get_ptr(m->sp_l)+2, get_ptr(m->iw));
    // Get rows in L (strictly lower entries only)
    int nnz_l = m->sp_l[2+n];
    m->sp_l.resize(2 + n+1 + nnz_l);
    casadi_ldl_row(sp_, get_ptr(m->parent), get_ptr(m->sp_l)+2, get_ptr(m->sp_l)+2+n+1,
                    get_ptr(m->iw));
    m->nz_l.resize(nnz_l);
    m->d.resize(n);
    m->w.resize(n);
    return 0;
  }

  int LinsolLdl::sfact(void* mem, const double* A) const {
    return 0;
  }

  int LinsolLdl::nfact(void* mem, const double* A) const {
    auto m = static_cast<LinsolLdlMemory*>(mem);
    casadi_ldl(sp_, get_ptr(m->parent), get_ptr(m->sp_l),
               A, get_ptr(m->nz_l), get_ptr(m->d), get_ptr(m->iw), get_ptr(m->w));
    return 0;
  }

  int LinsolLdl::solve(void* mem, const double* A, double* x, int nrhs, bool tr) const {
    auto m = static_cast<LinsolLdlMemory*>(mem);
    int n = this->nrow();
    for (int k=0; k<nrhs; ++k) {
      //      LDL'x = b <=> x = L\D\L'\b
      //  Solve for L'
      casadi_ldl_trs(get_ptr(m->sp_l), get_ptr(m->nz_l), x, 0);
      // Divide by D
      for (int i=0; i<n; ++i) x[i] /= m->d[i];
      // Solve for L
      casadi_ldl_trs(get_ptr(m->sp_l), get_ptr(m->nz_l), x, 1);
      // Next rhs
      x += n;
    }
    return 0;
  }

  int LinsolLdl::neig(void* mem) const {
    // Count number of negative eigenvalues
    auto m = static_cast<LinsolLdlMemory*>(mem);
    int n = this->nrow();
    int ret = 0;
    for (int i=0; i<n; ++i) if (m->d[i]<0) ret++;
    return ret;
  }

  int LinsolLdl::rank(void* mem) const {
    // Count number of nonzero eigenvalues
    auto m = static_cast<LinsolLdlMemory*>(mem);
    int n = this->nrow();
    int ret = 0;
    for (int i=0; i<n; ++i) if (m->d[i]!=0) ret++;
    return ret;
  }

} // namespace casadi
