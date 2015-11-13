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


#include "linsol.hpp"
#include "../std_vector_tools.hpp"
#include "../mx/mx_node.hpp"
#include <typeinfo>

using namespace std;
namespace casadi {

  Linsol::Linsol(const std::string& name, const Sparsity& sparsity, int nrhs)
    : FunctionInternal(name), sparsity_(sparsity), nrhs_(nrhs) {

    // Make sure arguments are consistent
    casadi_assert(!sparsity.isNull());
    casadi_assert_message(sparsity.size2()==sparsity.size1(),
                          "Linsol::init: the matrix must be square but got "
                          << sparsity.dim());
    casadi_assert_message(!sparsity.is_singular(),
                          "Linsol::init: singularity - the matrix is structurally "
                          "rank-deficient. sprank(J)=" << sprank(sparsity)
                          << " (in stead of "<< sparsity.size2() << ")");

    // Calculate the Dulmage-Mendelsohn decomposition
    std::vector<int> coarse_rowblock, coarse_colblock;
    sparsity.btf(rowperm_, colperm_, rowblock_, colblock_, coarse_rowblock, coarse_colblock);

    // Number of equations
    neq_ = sparsity.size2();

    ischeme_ = {"A", "B"};
    oscheme_ = {"X"};
  }

  Linsol::~Linsol() {
  }

  Sparsity Linsol::get_sparsity_in(int ind) const {
    switch (static_cast<LinsolInput>(ind)) {
    case LINSOL_A:
      return sparsity_;
    case LINSOL_B:
      return Sparsity::dense(neq_, nrhs_);
    case LINSOL_NUM_IN: break;
    }
    return Sparsity();
  }

  Sparsity Linsol::get_sparsity_out(int ind) const {
    switch (static_cast<LinsolOutput>(ind)) {
    case LINSOL_X:
      return Sparsity::dense(neq_, nrhs_);
    case LINSOL_NUM_OUT: break;
    }
    return Sparsity();
  }

  void Linsol::init() {
    // Call the base class initializer
    FunctionInternal::init();

  }

  void Linsol::eval(const double** arg, double** res, int* iw, double* w, void* mem) {
    // Get inputs and outputs
    const double *A = arg[LINSOL_A];
    const double *b = arg[LINSOL_B];
    double *x = res[LINSOL_X];

    // Create memory reference
    Memory m(this, arg, res, iw, w, mem);

    // A zero linear system would be singular
    casadi_assert(A!=0);

    // If output not requested, nothing to do
    if (!x) return;

    // If right hand side is zero, solution is trivially zero (if well-defined)
    if (!b) {
      casadi_fill(x, neq_*nrhs_, 0.);
      return;
    }

    // Factorize the linear system
    linsol_factorize(m, A);

    // Solve the factorized system
    casadi_copy(b, neq_*nrhs_, x);
    linsol_solve(m, x, nrhs_, false);
  }

  void Linsol::
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
      rhs[d] = tr ? B_hat - mul(A_hat.T(), X) : B_hat - mul(A_hat, X);
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

  void Linsol::
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

  void Linsol::
  linsol_spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, void* mem,
               bool tr, int nrhs) {
    // Sparsities
    const Sparsity& A_sp = sparsity_in(LINSOL_A);
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
      linsol_spsolve(X, tmp, tr);

      // Continue to the next right-hand-side
      B += n;
      X += n;
    }
  }

  void Linsol::
  linsol_spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, void* mem,
               bool tr, int nrhs) {
    // Sparsities
    const Sparsity& A_sp = sparsity_in(LINSOL_A);
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
      linsol_spsolve(tmp, X, !tr);

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

  void Linsol::
  linsol_spsolve(DM& X, const DM& B, bool tr) const {
    bvec_t* X_bvec = reinterpret_cast<bvec_t*>(X.ptr());
    const bvec_t* B_bvec = reinterpret_cast<const bvec_t*>(B.ptr());
    linsol_spsolve(X_bvec, B_bvec, tr);
  }

  void Linsol::linsol_spsolve(bvec_t* X, const bvec_t* B, bool tr) const {

    const Sparsity& A_sp = sparsity_in(LINSOL_A);
    const int* A_colind = A_sp.colind();
    const int* A_row = A_sp.row();
    int nb = rowblock_.size()-1; // number of blocks

    if (!tr) {
      for (int b=0; b<nb; ++b) { // loop over the blocks forward

        // Get dependencies from all right-hand-sides in the block ...
        bvec_t block_dep = 0;
        for (int el=rowblock_[b]; el<rowblock_[b+1]; ++el) {
          int rr = rowperm_[el];
          block_dep |= B[rr];
        }

        // ... as well as all other variables in the block
        for (int el=colblock_[b]; el<colblock_[b+1]; ++el) {
          int cc = colperm_[el];
          block_dep |= X[cc];
        }

        // Propagate ...
        for (int el=colblock_[b]; el<colblock_[b+1]; ++el) {
          int cc = colperm_[el];

          // ... to all variables in the block ...
          X[cc] |= block_dep;

          // ... as well as to other variables which depends on variables in the block
          for (int k=A_colind[cc]; k<A_colind[cc+1]; ++k) {
            int rr=A_row[k];
            X[rr] |= block_dep;
          }
        }
      }

    } else { // transpose
      for (int b=nb-1; b>=0; --b) { // loop over the blocks backward

        // Get dependencies ...
        bvec_t block_dep = 0;
        for (int el=colblock_[b]; el<colblock_[b+1]; ++el) {
          int cc = colperm_[el];

          // .. from all right-hand-sides in the block ...
          block_dep |= B[cc];

          // ... as well as from all depending variables ...
          for (int k=A_colind[cc]; k<A_colind[cc+1]; ++k) {
            int rr=A_row[k];
            block_dep |= X[rr];
          }
        }

        // Propagate to all variables in the block
        for (int el=rowblock_[b]; el<rowblock_[b+1]; ++el) {
          int rr = rowperm_[el];
          X[rr] |= block_dep;
        }
      }
    }
  }

  void Linsol::linsol_eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, void* mem,
                             bool tr, int nrhs) {
    casadi_error("Linsol::eval_sxLinsol not defined for class "
                 << typeid(*this).name());
  }

  MX Linsol::linsol_solve(const MX& A, const MX& B, bool tr) {
    return A->getSolve(B, tr, shared_from_this<Function>());
  }

  std::map<std::string, Linsol::Plugin> Linsol::solvers_;

  const std::string Linsol::infix_ = "linsol";

} // namespace casadi

