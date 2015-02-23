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


#include "linear_solver_internal.hpp"
#include "../std_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../mx/mx_node.hpp"
#include <typeinfo>

INPUTSCHEME(LinsolInput)
OUTPUTSCHEME(LinsolOutput)

using namespace std;
namespace casadi {

  LinearSolverInternal::LinearSolverInternal(const Sparsity& sparsity, int nrhs) {
    // Make sure arguments are consistent
    casadi_assert(!sparsity.isNull());
    casadi_assert_message(sparsity.size2()==sparsity.size1(),
                          "LinearSolverInternal::init: the matrix must be square but got "
                          << sparsity.dimString());
    casadi_assert_message(!sparsity.isSingular(),
                          "LinearSolverInternal::init: singularity - the matrix is structurally "
                          "rank-deficient. sprank(J)=" << sprank(sparsity)
                          << " (in stead of "<< sparsity.size2() << ")");

    // Calculate the Dulmage-Mendelsohn decomposition
    std::vector<int> coarse_rowblock, coarse_colblock;
    sparsity.dulmageMendelsohn(rowperm_, colperm_, rowblock_, colblock_,
                               coarse_rowblock, coarse_colblock);

    // Allocate inputs
    setNumInputs(LINSOL_NUM_IN);
    input(LINSOL_A) = DMatrix(sparsity);
    input(LINSOL_B) = DMatrix::zeros(sparsity.size2(), nrhs);

    // Allocate outputs
    setNumOutputs(LINSOL_NUM_OUT);
    output(LINSOL_X) = input(LINSOL_B);

    input_.scheme = SCHEME_LinsolInput;
    output_.scheme = SCHEME_LinsolOutput;
  }

  void LinearSolverInternal::init() {
    // Call the base class initializer
    FunctionInternal::init();

    // Not prepared
    prepared_ = false;
  }

  LinearSolverInternal::~LinearSolverInternal() {
  }

  void LinearSolverInternal::evaluate() {
    /*  Factorization fact;
        if (called_once) {
        // Check if any element has changed
        bool any_change = false;
        const vector<double>& val = input(0).data();
        for (int i=0; i<val.size(); ++i) {
        if (val[i] != a[i]) {
        any_change = true;
        break;
        }
        }

        // Reuse factored matrix if matrix hasn't changed
        fact = any_change ? SAMEPATTERN : FACTORED;
        } else {
        fact = DOFACT;
        called_once = true;
        }*/

    // Call the solve routine
    prepare();

    // Make sure preparation successful
    if (!prepared_)
      throw CasadiException("LinearSolverInternal::evaluate: Preparation failed");

    // Solve the factorized system
    solve(false);
  }

  void LinearSolverInternal::solve(bool transpose) {
    // Get input and output vector
    const vector<double>& b = input(LINSOL_B).data();
    vector<double>& x = output(LINSOL_X).data();
    int nrhs = input(LINSOL_B).size2();

    // Copy input to output
    copy(b.begin(), b.end(), x.begin());

    // Solve the factorized system in-place
    solve(getPtr(x), nrhs, transpose);
  }

  void LinearSolverInternal::evalFwdLinsol(const MX& X, const std::vector<cpv_MX>& fseed,
                                           const std::vector<pv_MX>& fsens, bool tr) {
    const MX& A = X->dep(1);
    std::vector<int> rhs_ind;
    std::vector<MX> rhs;
    std::vector<int> col_offset(1, 0);
    for (int d=0; d<fsens.size(); ++d) {
      const MX& B_hat = *fseed[d][0];
      const MX& A_hat = *fseed[d][1];

      // Get right hand side
      MX rhs_d;
      if (tr) {
        rhs_d = B_hat - mul(A_hat.T(), X);
      } else {
        rhs_d = B_hat - mul(A_hat, X);
      }

      // Simplifiy if zero
      if (rhs_d.isZero()) {
        *fsens[d][0] = MX(rhs_d.shape());
      } else {
        rhs.push_back(rhs_d);
        rhs_ind.push_back(d);
        col_offset.push_back(col_offset.back()+rhs_d.size2());
      }
    }

    if (!rhs.empty()) {
      // Solve for all directions at once
      rhs = horzsplit(solve(A, horzcat(rhs), tr), col_offset);

      // Save result
      for (int i=0; i<rhs.size(); ++i) {
        *fsens[rhs_ind[i]][0] = rhs[i];
      }
    }
  }

  void LinearSolverInternal::
  callAdjLinsol(const std::vector<MX>& arg, const std::vector<MX>& res,
                const std::vector<std::vector<MX> >& aseed,
                std::vector<std::vector<MX> >& asens, bool tr) {              
    // Number of derivatives
    int nadj = aseed.size();
    const MX& B = arg[0];
    const MX& A = arg[1];
    const MX& X = res[0];

    // Solve for all directions at once
    std::vector<MX> rhs(nadj);
    std::vector<int> col_offset(nadj+1, 0);
    for (int d=0; d<nadj; ++d) {
      rhs[d] = aseed[d][0];
      col_offset[d+1] = col_offset[d] + rhs[d].size2();
    }
    rhs = horzsplit(solve(A, horzcat(rhs), !tr), col_offset);

    // Collect sensitivities
    asens.resize(nadj);
    for (int d=0; d<nadj; ++d) {
      asens[d].resize(2);

      // Propagate to A
      if (!tr) {
        asens[d][1] = -mul(rhs[d], X.T(), MX::zeros(A.sparsity()));
      } else {
        asens[d][1] = -mul(X, rhs[d].T(), MX::zeros(A.sparsity()));
      }

      // Propagate to B
      asens[d][0] = rhs[d];
    }
  }

  void LinearSolverInternal::
  spFwdLinsol(const std::vector<const bvec_t*>& arg,
              const std::vector<bvec_t*>& res, int* itmp, bvec_t* rtmp,
              bool tr, int nrhs) {
    // Sparsities
    const Sparsity& A_sp = input(LINSOL_A).sparsity();
    const int* A_colind = A_sp.colind();
    const int* A_row = A_sp.row();
    int n = A_sp.size1();

    // Get pointers to data
    const bvec_t *B=arg[0], *A = arg[1];
    bvec_t* X = res[0];
    bvec_t* tmp = rtmp;

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
      spSolve(X, tmp, tr);

      // Continue to the next right-hand-side
      B += n;
      X += n;
    }
  }

  void LinearSolverInternal::
  spAdjLinsol(const std::vector<bvec_t*>& arg,
              const std::vector<bvec_t*>& res, int* itmp, bvec_t* rtmp,
              bool tr, int nrhs) {
    // Sparsities
    const Sparsity& A_sp = input(LINSOL_A).sparsity();
    const int* A_colind = A_sp.colind();
    const int* A_row = A_sp.row();
    int n = A_sp.size1();

    // Get pointers to data
    bvec_t *B=arg[0], *A=arg[1], *X=res[0];
    bvec_t* tmp = rtmp;

    // For all right-hand-sides
    for (int r=0; r<nrhs; ++r) {
      // Solve transposed
      std::fill(tmp, tmp+n, 0);
      spSolve(tmp, B, !tr);

      // Clear seeds
      std::fill(B, B+n, 0);

      // Propagate to X
      for (int i=0; i<n; ++i) X[i] |= tmp[i];

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

  void LinearSolverInternal::spSolve(DMatrix& X, const DMatrix& B, bool transpose) const {
    bvec_t* X_bvec = reinterpret_cast<bvec_t*>(X.ptr());
    const bvec_t* B_bvec = reinterpret_cast<const bvec_t*>(B.ptr());
    spSolve(X_bvec, B_bvec, transpose);
  }

  void LinearSolverInternal::spSolve(bvec_t* X, const bvec_t* B, bool transpose) const {

    const Sparsity& A_sp = input(LINSOL_A).sparsity();
    const int* A_colind = A_sp.colind();
    const int* A_row = A_sp.row();
    int nb = rowblock_.size()-1; // number of blocks

    if (!transpose) {
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

  void LinearSolverInternal::evalSXLinsol(const cpv_SXElement& arg, const pv_SXElement& res,
                                           int* itmp, SXElement* rtmp, bool tr, int nrhs) {
    casadi_error("LinearSolverInternal::evalSXLinsol not defined for class "
                 << typeid(*this).name());
  }

  MX LinearSolverInternal::solve(const MX& A, const MX& B, bool transpose) {
    return A->getSolve(B, transpose, shared_from_this<LinearSolver>());
  }

  void LinearSolverInternal::solve(double* x, int nrhs, bool transpose) {
    casadi_error("LinearSolverInternal::solve not defined for class " << typeid(*this).name());
  }

  std::map<std::string, LinearSolverInternal::Plugin> LinearSolverInternal::solvers_;

  const std::string LinearSolverInternal::infix_ = "linearsolver";

  void LinearSolverInternal::solveL(double* x, int nrhs, bool transpose) {
    casadi_error("LinearSolverInternal::solveL not defined for class "
                 << typeid(*this).name());
  }

  Sparsity LinearSolverInternal::getFactorizationSparsity(bool transpose) const {
    casadi_error("LinearSolverInternal::getFactorizationSparsity not defined for class "
                 << typeid(*this).name());
    return Sparsity();
  }

  DMatrix LinearSolverInternal::getFactorization(bool transpose) const {
    casadi_error("LinearSolverInternal::getFactorization not defined for class "
                 << typeid(*this).name());
    return DMatrix();
  }

} // namespace casadi

