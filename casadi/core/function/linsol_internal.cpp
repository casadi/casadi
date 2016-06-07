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


#include "linsol_internal.hpp"

using namespace std;
namespace casadi {

  Function linsol_new(const std::string& name, const std::string& solver,
                  const Sparsity& sp, int nrhs, const Dict& opts) {
    Linsol F(name + "_linsol", solver, sp, opts);
    MX A = MX::sym("A", sp);
    MX b = MX::sym("b", sp.size2(), nrhs);
    MX x = F.solve(A, b);
    return Function(name, {A, b}, {x}, {"A", "B"}, {"X"});
  }

  LinsolInternal::LinsolInternal(const std::string& name, const Sparsity& sparsity)
    : FunctionInternal(name), sparsity_(sparsity) {

    // Make sure arguments are consistent
    casadi_assert(!sparsity.is_null());
    casadi_assert_message(sparsity.size2()==sparsity.size1(),
                          "LinsolInternal::init: the matrix must be square but got "
                          << sparsity.dim());
    casadi_assert_message(!sparsity.is_singular(),
                          "LinsolInternal::init: singularity - the matrix is structurally "
                          "rank-deficient. sprank(J)=" << sprank(sparsity)
                          << " (in stead of "<< sparsity.size2() << ")");
  }

  LinsolInternal::~LinsolInternal() {
  }

  void LinsolInternal::init(const Dict& opts) {
    // Call the base class initializer
    FunctionInternal::init(opts);

  }

  void LinsolInternal::init_memory(void* mem) const {
    auto m = static_cast<LinsolMemory*>(mem);

    // Set initial sparsity pattern
    if (!sparsity_.is_null()) {
      m->sparsity = sparsity_.compress();
    }
  }

  void LinsolInternal::linsol_eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem,
                             bool tr, int nrhs) {
    casadi_error("eval_sx not defined for " + type_name());
  }

  void LinsolInternal::linsol_solve(void* mem, double* x, int nrhs, bool tr) const {
    casadi_error("'linsol_solve' not defined for " + type_name());
  }

  void LinsolInternal::linsol_solveL(void* mem, double* x, int nrhs, bool tr) const {
    casadi_error("'linsol_solveL' not defined for " + type_name());
  }

  void LinsolInternal::linsol_factorize(void* mem, const double* A) const {
    casadi_error("'linsol_factorize' not defined for " + type_name());
  }

  Sparsity LinsolInternal::linsol_cholesky_sparsity(void* mem, bool tr) const {
    casadi_error("'linsol_cholesky_sparsity' not defined for " + type_name());
  }

  DM LinsolInternal::linsol_cholesky(void* mem, bool tr) const {
    casadi_error("'linsol_cholesky' not defined for " + type_name());
  }

  std::map<std::string, LinsolInternal::Plugin> LinsolInternal::solvers_;

  const std::string LinsolInternal::infix_ = "linsol";

} // namespace casadi
