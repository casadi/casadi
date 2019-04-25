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


#include "lsqr.hpp"

#ifdef WITH_DL
#include <cstdlib>
#endif // WITH_DL

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINSOL_LSQR_EXPORT
  casadi_register_linsol_lsqr(LinsolInternal::Plugin* plugin) {
    plugin->creator = Lsqr::creator;
    plugin->name = "lsqr";
    plugin->doc = Lsqr::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Lsqr::options_;
    plugin->deserialize = &Lsqr::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_LSQR_EXPORT casadi_load_linsol_lsqr() {
    LinsolInternal::registerPlugin(casadi_register_linsol_lsqr);
  }

  Lsqr::Lsqr(const std::string& name, const Sparsity& sp) :
    LinsolInternal(name, sp) {
  }

  Lsqr::~Lsqr() {
    clear_mem();
  }

  int Lsqr::init_mem(void* mem) const {
    if (LinsolInternal::init_mem(mem)) return 1;
    auto m = static_cast<LsqrMemory*>(mem);

    // Temporary storage
    m->w.resize(nrow()+4*ncol());
    m->A.resize(sp_.nnz());
    return 0;
  }

  int Lsqr::nfact(void* mem, const double* A) const {
    auto m = static_cast<LsqrMemory*>(mem);

    std::copy(A, A+m->A.size(), get_ptr(m->A));
    return 0;
  }

  void Lsqr::generate(CodeGenerator& g, const std::string& A, const std::string& x,
                          casadi_int nrhs, bool tr) const {
    // Codegen the integer vectors
    string sp = g.sparsity(sp_);

        // Place in block to avoid conflicts caused by local variables
    g << "{\n";
    g.comment("FIXME(@jaeandersson): Memory allocation can be avoided");
    g << "casadi_real w[" << nrow() + 4*ncol() << "];\n";

    // Solve
    g << g.lsqr_solve(A, x, nrhs, tr, sp, "w") << "\n";

    // End of block
    g << "}\n";

  }

  int Lsqr::solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const {
    auto m = static_cast<LsqrMemory*>(mem);
    return casadi_lsqr_solve(A, x, nrhs, tr, sp_, get_ptr(m->w));
  }

} // namespace casadi
