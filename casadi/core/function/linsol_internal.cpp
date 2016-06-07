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
#include "../std_vector_tools.hpp"
#include "../mx/mx_node.hpp"
#include <typeinfo>

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

  LinsolInternal::LinsolInternal(const std::string& name, const Sparsity& sparsity, int nrhs)
    : FunctionInternal(name), sparsity_(sparsity), nrhs_(nrhs) {

    // Make sure arguments are consistent
    casadi_assert(!sparsity.is_null());
    casadi_assert_message(sparsity.size2()==sparsity.size1(),
                          "LinsolInternal::init: the matrix must be square but got "
                          << sparsity.dim());
    casadi_assert_message(!sparsity.is_singular(),
                          "LinsolInternal::init: singularity - the matrix is structurally "
                          "rank-deficient. sprank(J)=" << sprank(sparsity)
                          << " (in stead of "<< sparsity.size2() << ")");

    // Calculate the Dulmage-Mendelsohn decomposition
    btf_ = sparsity.btf();

    // Number of equations
    neq_ = sparsity.size2();
  }

  LinsolInternal::~LinsolInternal() {
  }

  void LinsolInternal::init(const Dict& opts) {
    // Call the base class initializer
    FunctionInternal::init(opts);

  }

  void LinsolInternal::
  linsol_forward(const std::vector<MX>& arg, const std::vector<MX>& res,
                 const std::vector<std::vector<MX> >& fseed,
                 std::vector<std::vector<MX> >& fsens, bool tr) {
    // Number of derivatives
    int nfwd = fseed.size();
    const MX& A = arg[1];
    const MX& X = res[0];

    // Solve for all directions at once
    std::vector<MX> rhs(nfwd);
    std::vector<int> col_offset(nfwd+1, 0);
    for (int d=0; d<nfwd; ++d) {
      const MX& B_hat = fseed[d][0];
      const MX& A_hat = fseed[d][1];
      rhs[d] = tr ? B_hat - mtimes(A_hat.T(), X) : B_hat - mtimes(A_hat, X);
      col_offset[d+1] = col_offset[d] + rhs[d].size2();
    }
    rhs = horzsplit(linsol_solve(A, horzcat(rhs), tr), col_offset);

    // Fetch result
    fsens.resize(nfwd);
    for (int d=0; d<nfwd; ++d) {
      fsens[d].resize(1);
      fsens[d][0] = rhs[d];
    }
  }

  void LinsolInternal::
  linsol_reverse(const std::vector<MX>& arg, const std::vector<MX>& res,
                 const std::vector<std::vector<MX> >& aseed,
                 std::vector<std::vector<MX> >& asens, bool tr) {
    // Number of derivatives
    int nadj = aseed.size();
    const MX& A = arg[1];
    const MX& X = res[0];

    // Solve for all directions at once
    std::vector<MX> rhs(nadj);
    std::vector<int> col_offset(nadj+1, 0);
    for (int d=0; d<nadj; ++d) {
      rhs[d] = aseed[d][0];
      col_offset[d+1] = col_offset[d] + rhs[d].size2();
    }
    rhs = horzsplit(linsol_solve(A, horzcat(rhs), !tr), col_offset);

    // Collect sensitivities
    asens.resize(nadj);
    for (int d=0; d<nadj; ++d) {
      asens[d].resize(2);

      // Propagate to A
      MX a;
      if (!tr) {
        a = -mac(rhs[d], X.T(), MX::zeros(A.sparsity()));
      } else {
        a = -mac(X, rhs[d].T(), MX::zeros(A.sparsity()));
      }
      if (asens[d][1].is_empty(true)) {
        asens[d][1] = a;
      } else {
        asens[d][1] += a;
      }

      // Propagate to B
      if (asens[d][0].is_empty(true)) {
        asens[d][0] = rhs[d];
      } else {
        asens[d][0] += rhs[d];
      }
    }
  }

  void LinsolInternal::
  linsol_sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem,
               bool tr, int nrhs) {
    // Sparsities
    const Sparsity& A_sp = sparsity_;
    const int* A_colind = A_sp.colind();
    const int* A_row = A_sp.row();
    int n = A_sp.size1();

    // Get pointers to data
    const bvec_t *B=arg[0], *A = arg[1];
    bvec_t* X = res[0];
    bvec_t* tmp = w;

    // For all right-hand-sides
    for (int r=0; r<nrhs; ++r) {
      // Copy B to a temporary vector
      copy(B, B+n, tmp);

      // Add A_hat contribution to tmp
      for (int cc=0; cc<n; ++cc) {
        for (int k=A_colind[cc]; k<A_colind[cc+1]; ++k) {
          int rr = A_row[k];
          tmp[tr ? cc : rr] |= A[k];
        }
      }

      // Propagate to X
      std::fill(X, X+n, 0);
      A_sp.spsolve(btf_, X, tmp, tr);

      // Continue to the next right-hand-side
      B += n;
      X += n;
    }
  }

  void LinsolInternal::
  linsol_sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem,
               bool tr, int nrhs) {
    // Sparsities
    const Sparsity& A_sp = sparsity_;
    const int* A_colind = A_sp.colind();
    const int* A_row = A_sp.row();
    int n = A_sp.size1();

    // Get pointers to data
    bvec_t *B=arg[0], *A=arg[1], *X=res[0];
    bvec_t* tmp = w;

    // For all right-hand-sides
    for (int r=0; r<nrhs; ++r) {
      // Solve transposed
      std::fill(tmp, tmp+n, 0);
      A_sp.spsolve(btf_, tmp, X, !tr);

      // Clear seeds
      std::fill(X, X+n, 0);

      // Propagate to B
      for (int i=0; i<n; ++i) B[i] |= tmp[i];

      // Propagate to A
      for (int cc=0; cc<n; ++cc) {
        for (int k=A_colind[cc]; k<A_colind[cc+1]; ++k) {
          int rr = A_row[k];
          A[k] |= tmp[tr ? cc : rr];
        }
      }

      // Continue to the next right-hand-side
      B += n;
      X += n;
    }
  }

  void LinsolInternal::linsol_eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem,
                             bool tr, int nrhs) {
    casadi_error("eval_sx not defined for " + type_name());
  }

  MX LinsolInternal::linsol_solve(const MX& A, const MX& B, bool tr) {
    return A->getSolve(B, tr, shared_from_this<Linsol>());
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
