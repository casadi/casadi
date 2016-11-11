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
    Linsol F(name + "_linsol", solver, opts);
    MX A = MX::sym("A", sp);
    MX b = MX::sym("b", sp.size2(), nrhs);
    MX x = F.solve(A, b);
    return Function(name, {A, b}, {x}, {"A", "B"}, {"X"});
  }

  LinsolInternal::LinsolInternal(const std::string& name) : FunctionInternal(name) {
  }

  LinsolInternal::~LinsolInternal() {
  }

  void LinsolInternal::init(const Dict& opts) {
    // Call the base class initializer
    FunctionInternal::init(opts);

  }

  void LinsolInternal::init_memory(void* mem) const {
    //auto m = static_cast<LinsolMemory*>(mem);
  }

  void LinsolInternal::linsol_eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem,
                             bool tr, int nrhs) {
    casadi_error("eval_sx not defined for " + type_name());
  }

  void LinsolInternal::solve(void* mem, double* x, int nrhs, bool tr) const {
    casadi_error("'solve' not defined for " + type_name());
  }

  void LinsolInternal::solve_cholesky(void* mem, double* x, int nrhs, bool tr) const {
    casadi_error("'solve_cholesky' not defined for " + type_name());
  }

  void LinsolInternal::reset(void* mem, const int* sp) const {
    auto m = static_cast<LinsolMemory*>(mem);

    // Decompress pattern
    int nrow = *sp++;
    int ncol = *sp++;
    const int* colind = sp;
    int nnz = colind[ncol];
    const int* row = nnz==nrow*ncol ? 0 : sp+ncol+1;

    // Save to sparsity field
    m->sparsity.clear();
    m->sparsity.push_back(nrow);
    m->sparsity.push_back(ncol);
    m->sparsity.insert(m->sparsity.end(), colind, colind+ncol+1);
    if (row) {
      m->sparsity.insert(m->sparsity.end(), row, row+nnz);
    } else {
      for (int cc=0; cc<ncol; ++cc) {
        for (int rr=0; rr<nrow; ++rr) {
          m->sparsity.push_back(rr);
        }
      }
    }
  }

  void LinsolInternal::factorize(void* mem, const double* A) const {
    casadi_error("'factorize' not defined for " + type_name());
  }

  Sparsity LinsolInternal::linsol_cholesky_sparsity(void* mem, bool tr) const {
    casadi_error("'linsol_cholesky_sparsity' not defined for " + type_name());
  }

  DM LinsolInternal::linsol_cholesky(void* mem, bool tr) const {
    casadi_error("'linsol_cholesky' not defined for " + type_name());
  }

  int LinsolInternal::neig(void* mem) const {
    casadi_error("'neig' not defined for " + type_name());
  }

  int LinsolInternal::rank(void* mem) const {
    casadi_error("'rank' not defined for " + type_name());
  }

  std::map<std::string, LinsolInternal::Plugin> LinsolInternal::solvers_;

  const std::string LinsolInternal::infix_ = "linsol";

} // namespace casadi
