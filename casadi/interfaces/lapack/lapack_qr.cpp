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
    plugin->version = CASADI_VERSION;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_LAPACKQR_EXPORT casadi_load_linsol_lapackqr() {
    LinsolInternal::registerPlugin(casadi_register_linsol_lapackqr);
  }

  LapackQr::LapackQr(const std::string& name) :
    LinsolInternal(name) {
  }

  LapackQr::~LapackQr() {
    clear_memory();
  }

  Options LapackQr::options_
  = {{&FunctionInternal::options_},
     {{"max_nrhs",
       {OT_INT,
        "Maximum number of right-hand-sides that get processed in a single pass [default:10]."}}
     }
  };

  void LapackQr::init(const Dict& opts) {
    // Call the base class initializer
    LinsolInternal::init(opts);

    max_nrhs_ = 10;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="max_nrhs") {
        max_nrhs_ = op.second;
      }
    }
  }

  void LapackQr::init_memory(void* mem) const {
    LinsolInternal::init_memory(mem);
  }

  void LapackQr::reset(void* mem, const int* sp) const {
    LinsolInternal::reset(mem, sp);
    auto m = static_cast<LapackQrMemory*>(mem);
    m->mat.resize(m->ncol() * m->ncol());
    m->tau.resize(m->ncol());
    m->work.resize(max(max_nrhs_, m->ncol())*10);
  }

  void LapackQr::factorize(void* mem, const double* A) const {
    auto m = static_cast<LapackQrMemory*>(mem);

    int ret = casadi_lapackqr_factorize(A, m->ncol(), m->work.size(), get_ptr(m->sparsity),
      get_ptr(m->mat), get_ptr(m->tau), get_ptr(m->work));
    casadi_assert_message(ret == 0, "LapackQr::prepare: dgeqrf_ "
                                      "failed to factorize the Jacobian. Info: " << ret << ".");

  }

  void LapackQr::solve(void* mem, double* x, int nrhs, bool tr) const {
    auto m = static_cast<LapackQrMemory*>(mem);

    int ret = casadi_lapackqr_solve(max_nrhs_, m->nrow(), nrhs, tr, x, m->ncol(), m->tau.size(),
      m->work.size(), get_ptr(m->mat), get_ptr(m->tau), get_ptr(m->work));
    if (ret) casadi_assert_message(ret == 0, "LapackQr::solve: dormqr_ failed "
                                          "to solve the linear system. Info: " << ret << ".");
  }

} // namespace casadi
